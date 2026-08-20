# The decode-remodulate reference — scoping (2026-08-20)

The first-look verdict promoted this from v2 polish to the entry ticket:
single-channel PBR needs a reference that is not the received signal, and
the one clean reference available to a one-SDR rig is the waveform we can
re-synthesize from our own decoder's output. Scoping study ran over the
whole atsc3 tree; condensed here. (Full agent report in the 8/20 session;
every claim below carries file:line in that report.)

## Verdict: buildable in ~1000 LOC / 10 files / 3–5 days

Most of the hard parts already exist and are spec-gated:

- **FEC-corrected codewords exist at runtime and are simply discarded**
  (`m9_fast.py:1124`, `m10_core.py:269`). No LDPC re-encode needed — the
  BICM inverse is one gather each way (`lam_of_q`), then
  `m6_bicm.map_cells()`. A full forward TX chain even exists
  (`m6_bicm.PlpChain.encode`), exercised on air.
- **The cell-domain decode→re-encode→cancel idea already has a working
  prototype**: `e69_ldm.py:97-113`, complete with the corrupt-one-bit
  negative control. We are promoting it one layer, to the waveform.
- **Pilot lattice fully parameterized** (`spec_pilots.py`, gated against
  both A/322 editions); pilot PRBS + phases gated (`m3_spec.py:271`);
  **bootstrap is already bit-exact synthesizable** (`bootstrap.py:341`).
- **Missing entirely** (the actual build): a forward frequency
  interleaver (~8 LOC — the inverse exists), the pilot **amplitude** map
  (~45 LOC — every demod path passes amp=1.0 because it divides pilots
  into H; the boost tables are banked and Table 9.17 gives a free gate),
  and an OFDM IFFT/CP synthesizer (nothing in the tree synthesizes a
  data symbol).

**Recommended shape (option c, hybrid):** pilots + dummy cells + SBS
nulls from spec parameters; data cells re-mapped from the corrected
codeword where FEC converged, hard decisions where not; L1 cells from the
re-encoded L1 codeword; preamble synthesized (fully parameterized,
`m3_preamble.py`); bootstrap via the existing synthesizer at 9/8
resample. New `lab/m45_remod.py` (~550 LOC) + `gate_e94_remod.py`
(~250 LOC) + 9 small taps.

The gates that make it evidence: round-trip `demod(synth(X)) == X`;
LS-fit-one-H residual on air must drop ≥N dB (pilots alone ≫20 dB); the
e69 corrupt-one-bit negative control; samples-in == samples-out.

## The three risks that can eat the project

1. **Raw-timebase bookkeeping.** The decode chain anchors, truncates,
   de-rotates CFO, notches (discarding 128 samples on one path and 0 on
   the other, chosen at runtime), and fine-times ±20 samples per frame.
   The reference must be rebuilt **on the raw capture's clock**, CFO
   restored, un-notched, re-anchored per frame (real TX/RX clock drift is
   documented in `FineTrack`). Budget a day; gate with a delay-0 CAF peak
   in bin 0 ± ½ sample.
2. **LDM caps data-cell cancellation at ~5 dB on RF25.** The data cells
   carry two superposed layers and we have never decoded the enhanced one
   (e69: 0/24 enhanced blocks converged; core-only cancellation measured
   3.4 dB). A core-only reference cannot cancel what it cannot describe.
   **The pilots, preamble, and bootstrap are single-layer and cancel
   fully — and the pilot lattice is exactly what the CAF currently
   sees.** So build pilot-only and pilot+data references, A/B them, and
   report cancellation per class. RF33 (no LDM) is the friendlier first
   target for full-waveform cancellation.
3. **Pilot amplitude ratios.** A wrong A_SP:A_CP:A_edge injects a
   periodic error on precisely the lattice we're trying to cancel —
   a residual that impersonates the artifact. Mitigate twice: gate the
   table against printed Table 9.17, and independently MEASURE the
   per-class ratio from air (|Y/H| at each pilot class); require
   agreement before trusting either.

Guard fields before assuming: `L1B_papr_reduction` (tone reservation
would make some carriers unknowable), `L1B_first_sub_miso` (unmodeled
precoder). Both 0 in every banked L1 so far. Transmitter time-domain
windowing is unmodeled and is the classic -35 dB cancellation floor;
measure symbol-edge residuals before blaming anything else.

## Step 0, before any code

Re-run `m8_l1.py` on the radar capture itself
(`rf25_fox_g4_0810_2254.cs16`). The newest banked RF25 L1 reads
**8K FFT / GI 512 / SP12_4 / Cred 3** — which is what re-derived the
first-look artifact correctly (170.7 samples = 8192/(Dx·Dy=48);
198.53 Hz = symbol rate 794.1 / Dy=4). Every geometry number downstream
flows from that decode. Never trust a remembered config over the frame's
own L1.

## Sample-rate note (good news)

Captures at 6.912 MS/s hit the decoder's no-op resample path — the
reference and surveillance share one clock and one index space. Keep
radar captures at 6.912 and the whole resampling risk class vanishes
(the 8.0 MS/s-era captures do NOT qualify).

## Step 0 — DONE 8/20 17:40Z, from the radar capture's own L1

`m8_l1.py rf25_fox_g4_0810_2254.cs16 --rate 6912000 --frames 3` →
`lab/l1_radar_capture.json`. Signalled, not remembered:

| field | value |
|---|---|
| FFT / GI | **8192 / 512** |
| scattered pilot pattern | index 9 = **SP12_4** → Dx=12, Dy=4 |
| Cred | 3 |
| SBS | first + last |
| scattered_pilot_boost | 1 |
| preamble_num_symbols | 1 |
| **papr_reduction** | **0** ← no tone reservation; resynthesis is clean |
| **miso** | **0** ← no unmodeled precoder |
| core PLP0 | QPSK 9/15 |
| enhanced PLP1 | 256QAM-NUC 7/15, **ldm_injection_level 10 = 5.0 dB** |

**Both artifact observables reproduce from the signalled config:**
- symbol rate = 6.912e6/(8192+512) = 794.118 Hz; lattice repeats every
  Dy=4 symbols → **198.529 Hz** (observed 195.3–199.2, bin width 3.9)
- within-symbol pilot spacing Dx·Dy = 48 carriers → N/48 = **170.67
  samples = 7.402 km** (observed 7.37, bin width 43.4 m)

Mechanism settled. The two risk flags the scoping study said to guard —
PAPR tone reservation and MISO — are both **0**, so exact resynthesis is
not blocked. LDM is confirmed at 5.0 dB injection with the enhanced layer
in 256QAM-NUC, so the data-cell cancellation cap stands and the
pilot-only reference is the right first build.
