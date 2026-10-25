"""vsb_trellis.py -- trellis-decoded symbol decisions for the 8-VSB reference.

A slicer decides each 8-VSB symbol on its own. At 15 dB equalized SNR (the
Baltimore stations from here) that is wrong about once in ten symbols, and
a reference wrong 10 % of the time cannot cancel more than ~10 dB. The
broadcast is trellis coded, and a Viterbi decoder uses that.

The code was identified from our own decisions on a clean capture (RF31,
29.5 dB, 1.8 M symbols), not taken on faith: of every linear relation among
a symbol's three bits and the two previous symbols of the same encoder, one
holds on 100.00 % of symbols and every other one on 50.1 %:

    12 interleaved encoders; data symbol j (0..827, after the 4 sync
    symbols) of the q-th segment after a field sync belongs to encoder
    (j + 4 q) mod 12; with level index 4 Z2 + 2 Z1 + Z0 (-7 .. +7),
    each encoder obeys  Z0[n] = Z1[n-1] xor Z0[n-2].

Z1 is the coded input; Z0 is set by the state, so each branch selects a
subset of two levels 8 apart ({Z1 Z0} with Z2 = 0 or 1), and Z2 is decided
inside the subset, where the levels are 4x further apart than for the slicer.
Encoders are not clocked during field-sync segments, so each encoder's
sequence runs on across them.
"""
import numpy as np
import numba

LV = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float32)
SEG = 832
FIELD = 313


@numba.njit(cache=True)
def viterbi(u, lv):
    """8-state Viterbi on one encoder's soft symbols; returns level indices.
    state = (Z1[n-1], Z0[n-1], Z0[n-2]) packed as z1p<<2 | z0p<<1 | z0pp."""
    n = len(u)
    pm = np.zeros(8)
    prev = np.empty((n, 8), np.uint8)
    lvl = np.empty((n, 8), np.uint8)
    for t in range(n):
        new = np.full(8, 1e30)
        x = u[t]
        for s in range(8):
            z1p = s >> 2; z0p = (s >> 1) & 1; z0pp = s & 1
            z0 = z1p ^ z0pp
            for z1 in range(2):
                l0 = z1 * 2 + z0
                d0 = (x - lv[l0]) ** 2
                d1 = (x - lv[l0 + 4]) ** 2
                if d0 <= d1:
                    m = pm[s] + d0; L = l0
                else:
                    m = pm[s] + d1; L = l0 + 4
                ns = (z1 << 2) | (z0 << 1) | z0p
                if m < new[ns]:
                    new[ns] = m; prev[t, ns] = s; lvl[t, ns] = L
        mn = new.min()
        for s in range(8):
            pm[s] = new[s] - mn
    out = np.empty(n, np.uint8)
    s = int(np.argmin(pm))
    for t in range(n - 1, -1, -1):
        out[t] = lvl[t, s]
        s = prev[t, s]
    return out


def data_layout(n_sym, m0, fmask):
    """Indices of every data symbol in time order and its encoder number, or
    None if the block holds no field sync (q unknown)."""
    k0 = (m0 // 2) % SEG
    starts = np.arange(k0, n_sym - SEG + 1, SEG)
    fs = np.array([fmask[s + 10] for s in starts])
    if not fs.any():
        return None
    q = np.full(len(starts), -1)
    last = None
    for i, f in enumerate(fs):
        if f:
            last = i; q[i] = -1
        elif last is not None:
            q[i] = i - last
    first = int(np.argmax(fs))
    for i in range(first):                     # before the first field sync: count back
        q[i] = FIELD - (first - i)
    keep = (q >= 1) & (q <= FIELD - 1)
    j = np.arange(4, SEG)
    D = (starts[keep, None] + j[None, :]).ravel()
    E = ((j[None, :] - 4 + 4 * q[keep, None]) % 12).ravel()
    return D, E


def decide(u, a, m0, fmask):
    """Trellis-decoded levels in place of the slicer's (sync and field-sync
    symbols left as they are). Returns (a_new, fraction of data symbols changed)."""
    lay = data_layout(len(a), m0, fmask)
    if lay is None:
        return a, 0.0
    D, E = lay
    out = a.copy()
    for e in range(12):
        ii = D[E == e]
        out[ii] = LV[viterbi(u[ii].astype(np.float32), LV)]
    return out, float(np.mean(out[D] != a[D]))
