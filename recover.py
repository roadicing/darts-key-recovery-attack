#!/usr/bin/env python3
"""
recover.py

Description: Recovers the DARTS-128 secret key from the accumulated
  signature statistics (key-adaptive bias separation + tiered oracle
  search) and repacks it into the reference sk format. Before
  searching, the script estimates whether the signature count is
  sufficient; if it is clearly insufficient, it prints the estimated
  requirement instead of searching.

Input:  argv[1] = optional data directory (default: current directory);
                  reads <dir>/acc_a.bin, acc_b.bin, pk.bin (the signature
                  count is estimated from the accumulators themselves)
        --compare-sk = also compare the result against <dir>/sk.bin and
                  print the md5 of both (last 32 bytes excluded)

Output: writes <dir>/recovered_sk.bin (sk.bin format, seed K zeroed)
"""
import argparse
import hashlib
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from math import comb
from statistics import NormalDist

import numpy as np


N, L, M, TAU = 512, 2, 1024, 30
STATES = np.array([-1, 0, 1], dtype=np.int32)
_ND = NormalDist()

# per-signature per-coordinate spread of the statistic
# (scheme constant: depends on TAU and the rounding, not the key)
SIG1 = 28.0

# search schedule: per-run candidate budgets for the medium/deep tiers and
# the number of stage-1 models carried into the search
BUDGET1, BUDGET2, MAX_MODELS = 5_000_000, 100_000_000, 8


def Q(t):
    return _ND.cdf(-t)


def Qinv(p):
    return -_ND.inv_cdf(p)


# sufficiency diagnosis
def expected_errors(a, sigma, score):
    """
    expected stage-1 decision errors. The cross-validation score measures
    this key's mu-misfit variance inflation directly:
        score = 1 + Var(mu - mu_hat) / sigma^2   (from odd/even prediction)
    so the effective noise is sigma*sqrt(score). No empirical derating
    constant: the misfit that a global constant would approximate is measured
    per key, from the data. (score from half-fits slightly overestimates the
    full-fit misfit, which is the safe direction for a sufficiency test.)
    """
    if a <= 0 or sigma <= 0:
        return float('inf')
    sig_eff = sigma * np.sqrt(max(score, 1.0))
    return M * Q(0.5 * a / sig_eff)


def posterior_errors(mdl):
    """
    expected misclassified positions under the model's own posterior.
    """
    return float(np.sum(1 - mdl['prob'].max(axis=1)))


def search_capacity(n_suspects, budget):
    """
    errors guaranteed correctable within ONE dfs run over n_suspects with
    the given budget: the largest e with 1 + sum_{i<=e} C(S,i)*2^i <= budget
    (the 1 is the base candidate itself). Worst case over orderings; margin
    ordering does strictly better in practice.
    """
    total, e = 1, 0
    while e < n_suspects:
        cost = comb(n_suspects, e + 1) * (2 ** (e + 1))
        if total + cost > budget:
            break
        total += cost
        e += 1
    return e


def required_n(a, sigma, score, n_cur, target):
    """
    N at which the expected error count drops to `target`.
    """
    sig_eff = sigma * np.sqrt(max(score, 1.0))
    snr_cur = 0.5 * a / sig_eff
    if snr_cur <= 0:
        return float('inf')
    snr_req = Qinv(target / M)
    return n_cur * (snr_req / snr_cur) ** 2


def load_statistics(directory):
    acc = [np.fromfile(os.path.join(directory, nm),
                       dtype=np.int64).astype(np.float64)
           for nm in ('acc_a.bin', 'acc_b.bin')]
    if any(a.size != M for a in acc):
        raise ValueError(
            'each accumulator must contain exactly 1024 int64 values')
    # the signature count is estimated from the data itself: the odd/even
    # half-difference is pure noise with per-coordinate variance N*SIG1^2,
    # where SIG1 is the per-signature statistic spread (a scheme constant)
    d = (acc[0] - acc[1]) / TAU
    number = int(round(float(np.var(d, ddof=1)) / SIG1 ** 2))
    if number < 1:
        raise ValueError('accumulators look empty (no signatures collected)')
    # the halves hold the odd- and even-counter samples of the contiguous
    # counter range [0, N) that collect always sums completely, so their
    # sizes differ by at most one sample (odd N)
    na, nb = (number + 1) // 2, number // 2
    xa, xb = -acc[0] / (TAU * na), -acc[1] / (TAU * nb)
    v = float(np.var(xa - xb, ddof=1))
    if not np.isfinite(v) or v <= 0:
        raise ValueError(
            'split difference has no positive finite noise variance')
    sa = np.sqrt(v * nb / number)
    sb = np.sqrt(v * na / number)
    sig = np.sqrt(v * na * nb / number ** 2)
    x = (na * xa + nb * xb) / number
    z = float(np.mean(xa * xb)) * np.sqrt(M) / (sa * sb)
    return xa, xb, x, sa, sb, sig, z, number


# amplitude (a0) prior
def highpass(v, K):
    """
    remove the DC component and the K lowest Fourier modes per 512-block.
    """
    out = np.empty_like(v)
    for b in range(L):
        blk = v[b * N:(b + 1) * N]
        f = np.fft.rfft(blk)
        f[:K + 1] = 0
        out[b * N:(b + 1) * N] = np.fft.irfft(f, N)
    return out


def estimate_amplitude(xa, xb):
    """
    cross-spectrum plateau: E[hp(xa).hp(xb)] = a^2 * E[s^2] * frac_kept.
    """
    esp2 = 2 * 0.15  # E[s_j^2] = 2P with public P=0.15
    table = []
    for K in (2, 4, 8, 16, 24, 32, 48, 64, 96, 128):
        ha, hb = highpass(xa, K), highpass(xb, K)
        frac = (N - 2 * K - 1) / N if 2 * K + 1 < N else 0.0
        if frac <= 0:
            continue
        a2 = float(np.mean(ha * hb)) / (esp2 * frac)
        table.append((K, np.sqrt(max(a2, 0.0))))
    plateau = [a for K, a in table if K >= 16 and a > 0]
    a0 = float(np.median(plateau)) if plateau else 0.0
    return a0, table


# smoothers
def _med_1d(v, w, pad='reflect'):
    h = w // 2
    if pad == 'reflect':
        ext = np.concatenate([v[1:h + 1][::-1], v, v[-h - 1:-1][::-1]])
    else:  # circular
        ext = np.concatenate([v[-h:], v, v[:h]])
    return np.median(np.lib.stride_tricks.sliding_window_view(ext, w), axis=-1)


_BASIS = {}


def _lowk_1d(K):
    if K not in _BASIS:
        t = np.arange(N)
        cols = [np.ones(N) / np.sqrt(N)]
        for k in range(1, K + 1):
            cols.append(np.cos(2 * np.pi * k * t / N) * np.sqrt(2 / N))
            cols.append(np.sin(2 * np.pi * k * t / N) * np.sqrt(2 / N))
        _BASIS[K] = np.column_stack(cols)
    B = _BASIS[K]
    return lambda r: B @ (B.T @ r)


def _combo_1d(K, w):
    base = _lowk_1d(K)
    return lambda r: (lambda b: b + _med_1d(r - b, w, 'reflect'))(base(r))


def smoother_pool():
    """
    per-block smoothers: (tag, fn on 512-vector). Both reflect- and
    circular-padded medians are included on purpose: mu is not circular, but
    padding-mode diversity makes the models disagree exactly where the mu
    estimate is contaminated by edge effects or jumps.
    """
    pool = [(f'med{w}', (lambda w=w: lambda r: _med_1d(r, w, 'reflect'))())
            for w in (3, 5, 7, 9, 11, 15, 21, 31, 45, 63)]
    pool += [(f'cmed{w}', (lambda w=w: lambda r: _med_1d(r, w, 'circular'))())
             for w in (9, 31, 63)]
    pool += [(f'lowK{K}', _lowk_1d(K)) for K in (4, 8, 16, 32)]
    pool += [('k8m9', _combo_1d(8, 9)), ('k16m25', _combo_1d(16, 25))]
    return pool


def per_block(fn_pair):
    f0, f1 = fn_pair
    return lambda r: np.concatenate([f0(r[:N]), f1(r[N:])])


# soft-EM (GMM)
def posterior(resid, a, rho, sigma):
    z = np.stack([resid + a, resid, resid - a], axis=1)
    ll = -0.5 * (z / sigma) ** 2
    ll += np.log([rho / 2, 1 - rho, rho / 2])
    ll -= ll.max(axis=1, keepdims=True)
    p = np.exp(ll)
    return p / p.sum(axis=1, keepdims=True)


def fit_model(values, sigma, smooth, a0=None, max_iter=100):
    """
    damped soft-EM over the ternary mixture; smoother supplied as fn.
    """
    mu = smooth(values)
    resid = values - mu
    if a0 is not None and np.isfinite(a0) and a0 > 2 * sigma:
        a = float(a0)
    else:
        a = max(float(np.quantile(np.abs(resid), 0.9)), 2 * sigma)
    floor = 1 / (len(values) + 2)
    rho = float(np.clip(np.mean(np.abs(resid) > a / 2), floor, 1 - floor))
    converged, it = False, 0
    for it in range(1, max_iter + 1):
        p = posterior(values - mu, a, rho, sigma)
        ms = p[:, 2] - p[:, 0]
        ms2 = p[:, 2] + p[:, 0]
        new_a = max(float(np.dot(values - mu, ms) / max(ms2.sum(), 1e-12)),
                    0.1 * sigma)
        new_rho = float(np.clip(ms2.mean(), floor, 1 - floor))
        new_mu = smooth(values - new_a * ms)
        change = max(float(np.max(np.abs(new_mu - mu))) / sigma,
                     abs(new_a - a) / sigma, abs(new_rho - rho))
        mu = (mu + new_mu) / 2  # damping against median-update oscillation
        a, rho = new_a, new_rho
        if change < 1e-5:
            converged = True
            break
    p = posterior(values - mu, a, rho, sigma)
    s = STATES[np.argmax(p, axis=1)]
    qs = np.sort(p, axis=1)
    return {'s': s, 'mu': mu, 'a': a, 'rho': rho, 'prob': p,
            'margin': qs[:, -1] - qs[:, -2],
            'mean': mu + a * (p[:, 2] - p[:, 0]),
            'converged': converged, 'iters': it}


def select_models(xa, xb, x, sa, sb, sig, a0, rel_thresh=1.07, max_pairs=36):
    """
    cross-validate each smoother per block on the odd/even split, then
    fully fit the best per-block combinations plus every single smoother.
    """
    pool = smoother_pool()
    best_per_block = []
    all_scores = []  # per block: full {tag: score}, used for honest ranking
    for b in range(L):
        sl = slice(b * N, (b + 1) * N)
        scores = []
        for tag, fn in pool:
            fa = fit_model(xa[sl], sa, fn, a0)
            fb = fit_model(xb[sl], sb, fn, a0)
            sc = 0.5 * (np.mean((fa['mean'] - xb[sl]) ** 2) / sb ** 2 +
                        np.mean((fb['mean'] - xa[sl]) ** 2) / sa ** 2)
            scores.append((sc, tag, fn))
        scores.sort(key=lambda t: t[0])
        all_scores.append({tag: sc for sc, tag, _ in scores})
        thresh = scores[0][0] * rel_thresh
        keep = [s for s in scores if s[0] <= thresh][:6]
        if len(keep) < 3:
            keep = scores[:3]
        best_per_block.append(keep)

    models, seen = [], set()

    def add(mdl, tag, score):
        key = mdl['s'].tobytes()
        if key in seen:
            return
        seen.add(key)
        mdl['tag'] = tag
        mdl['score'] = score
        models.append(mdl)

    pairs = list(itertools.product(*best_per_block))
    pairs.sort(key=lambda pr: pr[0][0] + pr[1][0])
    for (s0, t0, f0), (s1, t1, f1) in pairs[:max_pairs]:
        add(fit_model(x, sig, per_block((f0, f1)), a0), f'{t0},{t1}',
            (s0 + s1) / 2)
    # single smoothers applied uniformly: structural diversity for the
    # suspect union and the disagreement vote; ranked by their real
    # cross-validation scores (no placeholder values)
    for tag, fn in dict(pool).items():
        sc = float(np.mean([all_scores[b][tag] for b in range(L)]))
        add(fit_model(x, sig, per_block((fn, fn)), a0), f'{tag},{tag}', sc)

    models.sort(key=lambda m: m['score'])
    return models


def sufficiency(a, sigma, score, n_cur, post_est, target, e_cap):
    """
    expected stage-1 errors vs search capacity. Two imperfect
    indicators are combined: the analytic model can miss mu jumps, the
    posterior is overconfident on confidently-wrong positions. Use the
    pessimistic one for the go-ahead and the optimistic one for giving up:
    SUFFICIENT only when both agree it is safe, INSUFFICIENT only when even
    the optimistic indicator says the search is hopeless.
    """
    e_form = expected_errors(a, sigma, score)
    e_opt, e_pes = min(e_form, post_est), max(e_form, post_est)
    n_req = required_n(a, sigma, score, n_cur, target)
    if e_pes <= max(1.0, e_cap / 2):
        verdict = 'SUFFICIENT'
    elif e_opt > 2.0 * e_cap:
        verdict = 'INSUFFICIENT'
    else:
        verdict = 'BORDERLINE'
    return e_form, e_opt, e_pes, n_req, verdict


# oracle DFS
def run_oracle(D, dfs_exe, base, suspects, maxdepth, budget):
    """
    one dfs run over the suspect set; returns (rec, ee, tried) with
    rec/ee None when the search limits (depth or budget) were reached
    without a hit.
    """
    dfs_exe = os.path.abspath(dfs_exe)
    if not os.path.isfile(dfs_exe):
        raise FileNotFoundError(
            f'dfs binary not found: {dfs_exe} (run build.sh first)')
    tmp = tempfile.mkdtemp(prefix='recover_')
    try:
        bp = os.path.join(tmp, 'base.bin')
        sp = os.path.join(tmp, 'susp.bin')
        np.asarray(base, dtype=np.int32).tofile(bp)
        np.asarray(suspects, dtype=np.int32).tofile(sp)
        r = subprocess.run([dfs_exe,
                            os.path.abspath(os.path.join(D, 'pk.bin')),
                            bp, sp, str(maxdepth), str(budget)],
                           capture_output=True, text=True, cwd=tmp)
        mt = re.search(r'tried=(\d+)', r.stdout)
        tried = int(mt.group(1)) if mt else 0
        if 'FOUND depth' in r.stdout:
            rec = np.fromfile(os.path.join(tmp, 'dfs_found.bin'),
                              dtype=np.int32)
            ee = np.fromfile(os.path.join(tmp, 'dfs_found_e.bin'),
                             dtype=np.int32)
            if (rec.size != M or ee.size != N
                    or not np.isin(rec, (-1, 0, 1)).all()):
                raise RuntimeError('dfs returned malformed result files')
            return rec, ee, tried
        if r.returncode not in (0, 2):  # dfs: 0 = found, 2 = limits reached
            err = r.stderr.strip()[:300]
            raise RuntimeError(
                f'dfs exited with code {r.returncode}: {err}')
        return None, None, tried
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def poly_p_pack_py(coeffs):
    t = (1 - np.asarray(coeffs, dtype=np.int32)) & 0x3
    t = t.reshape(-1, 4)
    packed = t[:, 0] | (t[:, 1] << 2) | (t[:, 2] << 4) | (t[:, 3] << 6)
    return packed.astype(np.uint8).tobytes()


# main
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('data_dir', nargs='?', default='.',
                    help='directory with acc_a.bin/acc_b.bin/pk.bin '
                         '(default: current directory)')
    ap.add_argument('--compare-sk', action='store_true',
                    help='compare the recovered key against '
                         '<data_dir>/sk.bin and print the md5 of both')
    args = ap.parse_args()
    D = args.data_dir
    dfs_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dfs')

    xa, xb, x, sa, sb, sig, z, Nsig = load_statistics(D)
    print(f'signatures N ~ {Nsig:,}  (estimated from the accumulators)')
    print(f'leakage detection: z = {z:+.2f}  (z > 3 means significant)')
    if z <= 3.0:
        print('signal not significant; collect more signatures.')
        return 2

    print('stage 1 bias separation ... ', end='', flush=True)
    a0, _ = estimate_amplitude(xa, xb)
    models = select_models(xa, xb, x, sa, sb, sig, a0)
    print('done', flush=True)

    a_est = float(np.median([m['a'] for m in models[:3]]))
    post_est = float(np.median([posterior_errors(m) for m in models[:3]]))
    # mu-misfit inflation of THIS key, measured by cross-validation
    score_est = float(np.median([m['score'] for m in models[:3]]))
    # correctable-error guarantee, computed honestly per round: each dfs call
    # re-enumerates from scratch, so rounds do not compose; the guarantee is
    # the best single-round full-coverage depth over that round's own
    # suspect-set size and budget. (round tiers defined below)
    tier_specs = [(96, BUDGET1), (96, BUDGET2), (192, 2 * BUDGET2)]
    e_cap = max(search_capacity(s, b) for s, b in tier_specs)
    target = max(1.0, e_cap / 2)
    e_form, e_opt, e_pes, n_req, verdict = sufficiency(
        a_est, sig, score_est, Nsig, post_est, target, e_cap)
    if verdict == 'INSUFFICIENT':
        if np.isfinite(n_req):
            print(f'signature count insufficient: estimated need '
                  f'~{n_req:,.0f} (collect ~{max(n_req - Nsig, 0):,.0f} '
                  f'more: rerun ./collect with the larger count and '
                  f'--load ckpt.bin)')
        else:
            print('signature count insufficient: amplitude too weak to '
                  'extrapolate a requirement; collect much more data')
        return 2

    # suspect set
    # suspicion(j) = max( mean over models of posterior doubt, fraction of
    # models disagreeing with the plurality vote ). The disagreement term is
    # essential: a mu jump makes every smoother *confidently* wrong (margin
    # 1.0), but different smoother families are wrong at *different* places,
    # so the vote exposes what the margins hide. A small per-block edge
    # allowance is added because median windows are least reliable there.
    # Disagreement/edge positions are placed first, so the cap evicts
    # low-suspicion positions before them (if the reserved set itself
    # exceeds the cap, its least-suspicious members are still dropped).
    def build_suspects(models, cap):
        ref = models[0]['s']
        votes = []
        for m in models:
            si = m['s']
            if np.count_nonzero(-si - ref) < np.count_nonzero(si - ref):
                si = -si
            votes.append(si)
        votes = np.asarray(votes)
        n_pos = (votes == 1).sum(0)
        n_neg = (votes == -1).sum(0)
        n_zer = (votes == 0).sum(0)
        maj = np.where(n_zer >= np.maximum(n_pos, n_neg), 0,
                       np.where(n_pos >= n_neg, 1, -1))
        disagree = (votes != maj).mean(axis=0)
        doubt = np.mean([1 - m['margin'] for m in models], axis=0)
        suspicion = np.maximum(doubt, disagree)
        reserved = set(np.flatnonzero(disagree > 0).tolist())
        for b in range(L):
            reserved |= set(range(b * N, b * N + 8))
            reserved |= set(range((b + 1) * N - 8, (b + 1) * N))
        cand = set(reserved)
        for m in models:
            cand |= set(np.argsort(m['margin'])[:40].tolist())
        cand |= set(np.argsort(-suspicion)[:cap].tolist())
        order = sorted(cand, key=lambda j: -suspicion[j])
        keep = [j for j in order if j in reserved][:cap]
        keep += [j for j in order if j not in reserved][:cap - len(keep)]
        return np.array(sorted(keep, key=lambda j: -suspicion[j]),
                        dtype=np.int32)

    suspects = build_suspects(models, 96)

    search_models = models[:MAX_MODELS]
    rec = ee = None
    # Shallow tiers rotate across ALL models before any model gets a deep
    # budget: the score ranking is imperfect, and on test_data the 2-error
    # model ranks 2nd while the 5-error 1st-ranked model would otherwise
    # burn 10M candidates at depths it cannot complete. Depth 2 over 96
    # suspects costs ~18k candidates per model-sign, so the sweep is cheap.
    # tier list: (label, models, suspect_cap, maxdepth, budget per run)
    tiers = [
        ('tier 1', search_models, 96, 2, 250_000),
        ('tier 2', search_models, 96, 5, BUDGET1),
        ('tier 3', search_models[:3], 96, 10, BUDGET2),
        ('tier 4', search_models[:2], 192, 14, 2 * BUDGET2),
    ]
    print('oracle search ... ', end='', flush=True)
    tried_total, hit_depth = 0, None
    for label, tier_models, cap, depth, budget in tiers:
        if rec is not None:
            break
        tier_suspects = suspects if cap == 96 else build_suspects(models, cap)
        for m in tier_models:
            for base in (m['s'], -m['s']):
                r_, e_, tried = run_oracle(D, dfs_exe, base, tier_suspects,
                                           depth, budget)
                tried_total += tried
                if r_ is not None:
                    rec, ee, hit_depth = r_, e_, depth
                    break
            if rec is not None:
                break
        if rec is None:
            print(f'not found at depth<={depth}; ', end='', flush=True)
    if rec is not None:
        print(f'key found after {tried_total:,} candidates '
              f'(depth<={hit_depth})')
    else:
        print('key not found within the enumeration budget')
        if np.isfinite(n_req) and n_req > Nsig:
            print(f'signature count insufficient: estimated need '
                  f'~{n_req:,.0f} (collect ~{n_req - Nsig:,.0f} more: '
                  f'rerun ./collect with the larger count and '
                  f'--load ckpt.bin)')
        else:
            print('diagnosis: the estimate predicted few errors yet the '
                  'search failed; doubling the signature count is the '
                  'cheap fix.')
        return 2

    pk = open(os.path.join(D, 'pk.bin'), 'rb').read()
    rebuilt = pk + poly_p_pack_py(rec[:N]) + poly_p_pack_py(rec[N:]) + \
        poly_p_pack_py(ee[-N:]) + bytes(32)
    out = os.path.normpath(os.path.join(D, 'recovered_sk.bin'))
    open(out, 'wb').write(rebuilt)
    print(f'recovered key written to {out}')

    if args.compare_sk:
        sk_path = os.path.normpath(os.path.join(D, 'sk.bin'))
        if not os.path.isfile(sk_path):
            print(f'note: no sk.bin found in {os.path.abspath(D)}; '
                  f'saved without comparison.')
        else:
            sk = open(sk_path, 'rb').read()
            if len(sk) != len(rebuilt):
                print(f'warning: {sk_path} has length {len(sk)}, expected '
                      f'{len(rebuilt)}; cannot compare')
                return 1
            cut = len(sk) - 32
            m1 = hashlib.md5(sk[:cut]).hexdigest()
            m2 = hashlib.md5(rebuilt[:cut]).hexdigest()
            print(f'comparing with {sk_path} (md5 of bytes [0:{cut}]; '
                  f'last 32 bytes = signing seed K, excluded):')
            print(f'  sk.bin           = {m1}\n  recovered_sk.bin = {m2}  '
                  f'{"<-- identical" if m1 == m2 else "<-- MISMATCH!"}')
            if m1 != m2:
                return 1
            print(f'\nPOC PASSED: the DARTS-128 secret key was fully '
                  f'recovered from ~{Nsig:,} public signatures alone.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
