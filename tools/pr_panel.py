"""pr_panel.py -- the range-Doppler scope. v1.0 2026-08-21, port 8650.

What it shows, per dwell: the flattened, CFAR-normalised range-Doppler map
as a heatmap; every CFAR detection as a ring; every ADS-B aircraft that
geometry puts inside the map as a labelled crosshair AT WHERE ITS OWN
TELEMETRY SAYS IT MUST APPEAR; the pilot comb teeth as faint lines; and
the kinematic tracks as polylines. A slider (and play) walks the dwells,
so a real target is a ring riding a crosshair, frame after frame.

That overlay is the instrument's whole honesty story in one picture:
the crosshairs are drawn from truth the radar never saw, so a ring on a
crosshair is evidence and a ring off every crosshair is residue.

House rules: localhost only; the observer position is read server-side
from the private GPS fix and NEVER sent to the browser - only predicted
(range, Doppler) pairs go over the wire; append-only, deletes nothing;
no clickable DOM is built from a render loop (8/04 law) - the slider and
run selector are static, the canvas is repainted.

Data it reads (all under lab/, produced by track_hunt.py --save-maps and
truth_capture.py):
    tracks_<run>.json          detections + tracks (track_hunt)
    maps/tracks_<run>/rd_*.npy per-dwell SNR maps, axes.npy
    truth_rf25_<stamp>.meta.json + .adsb.jsonl  (truth_capture)

Run:  python tools/pr_panel.py [--port 8650]
Open: http://127.0.0.1:8650/
"""
import argparse
import glob
import json
import os
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.join(os.path.dirname(HERE), "lab")
sys.path.insert(0, HERE)
import truth_compare as TC                                    # noqa: E402
import pr_pilots as PP                                        # noqa: E402

BUILD = "pr_panel 1.0 2026-08-21"
_cache = {}
_lock = threading.Lock()


# ----------------------------------------------------------------- data
def runs():
    out = []
    for p in sorted(glob.glob(os.path.join(LAB, "tracks_*.json"))):
        stem = os.path.splitext(os.path.basename(p))[0]
        maps = os.path.join(LAB, "maps", stem)
        out.append(dict(name=stem, tracks=p,
                        has_maps=os.path.isdir(maps)
                        and bool(glob.glob(os.path.join(maps, "rd_*.npy")))))
    return out


def find_meta(capture_path):
    """truth_capture's meta for a capture, if this capture came from it."""
    stem = os.path.splitext(os.path.basename(capture_path))[0]
    p = os.path.join(LAB, stem + ".meta.json")
    return json.load(open(p)) if os.path.exists(p) else None


def load_run(name):
    trk = json.load(open(os.path.join(LAB, name + ".json")))
    det = trk["detections"]
    times = sorted({d["t_s"] for d in det})
    maps_dir = os.path.join(LAB, "maps", name)
    map_files = {}
    if os.path.isdir(maps_dir):
        for f in glob.glob(os.path.join(maps_dir, "rd_*.npy")):
            t = float(os.path.basename(f)[3:-4])
            map_files[round(t, 4)] = f
        for t in sorted(map_files):
            if round(t, 4) not in {round(x, 4) for x in times}:
                times.append(t)
        times = sorted(set(round(t, 4) for t in times))
    axes = None
    ax_p = os.path.join(maps_dir, "axes.npy")
    if os.path.exists(ax_p):
        a = np.load(ax_p, allow_pickle=True)
        axes = dict(range_m=[float(x) for x in a[0]],
                    doppler_hz=[float(x) for x in a[2]])

    # ADS-B predictions, server-side only; observer never leaves this process
    preds = {}
    aircraft = []
    meta = find_meta(trk["capture"])
    if meta and os.path.exists(meta.get("truth", "")):
        truth = TC.load_truth(meta["truth"])
        lat, lon, alt = TC.observer()
        R = TC.ecef(lat, lon, alt)
        T = TC.ecef(TC.TX_LAT, TC.TX_LON, TC.TX_ANT_AMSL_M)
        dRT = float(np.linalg.norm(T - R))
        t0 = float(meta["t_launch_unix"])
        for icao, s in truth.items():
            cs = (s[0][7] or icao).strip()
            rows = []
            for ts in times:
                p = TC.predict(s, t0 + ts, T, R, dRT)
                if p and 0 < p[0] < 60e3 and abs(p[1]) < 1000:
                    rows.append([ts, round(p[0], 1), round(p[1], 2)])
            if rows:
                aircraft.append(dict(icao=icao, callsign=cs, n=len(rows)))
                preds[icao] = rows
    g_nfft, g_pat = 8192, "SP12_4"
    try:
        l1 = json.load(open(os.path.join(LAB, "l1_radar_capture.json")))
        g_nfft = int(l1["geometry"]["nfft"])
        g_pat = PP.P.PATTERNS[int(l1["l1_basic"]
                                  ["L1B_first_sub_scattered_pilot_pattern"])] \
            if hasattr(PP, "P") else g_pat
    except Exception:
        pass
    tooth = PP.lattice_range_m(g_nfft, g_pat)
    run = dict(name=name, build=BUILD, capture=os.path.basename(trk["capture"]),
               times=times, detections=det, tracks=trk.get("tracks", []),
               comb_m=[tooth * n for n in range(1, 9)],
               axes=axes, has_maps=bool(map_files),
               aircraft=sorted(aircraft, key=lambda a: -a["n"]), preds=preds,
               note=trk.get("note", ""))
    # No caching. A run loaded while track_hunt --save-maps was still
    # writing got frozen at a partial map set for the panel's whole life
    # (8/21, 7 of 54 maps). Rebuilding costs a glob and ~60 aircraft x 54
    # dwells of predict() - milliseconds - and is always current.
    return run, map_files


def map_json(name, t_s, ds_r=2, ds_d=1):
    run, map_files = load_run(name)
    key = min(map_files, key=lambda t: abs(t - t_s)) if map_files else None
    if key is None or abs(key - t_s) > 0.2:
        return None
    m = np.load(map_files[key]).astype(np.float32)
    m = m[::ds_r, ::ds_d]
    # The map is CFAR-normalised dB (signal over LOCAL noise). The floor is
    # ~0 dB by construction and must render DARK - a percentile stretch
    # painted the speckle bright amber and buried the exceedances (8/21,
    # first screenshot). Fixed stretch: floor -> black, +16 dB -> white,
    # gamma to keep the first few dB quiet.
    lo = float(np.median(m))
    hi = max(lo + 16.0, float(np.percentile(m, 99.98)))
    q = np.clip((m - lo) / max(hi - lo, 1e-6), 0, 1) ** 1.8
    q = (q * 255).astype(np.uint8)
    return dict(t_s=key, shape=list(q.shape), lo_db=round(lo, 1),
                hi_db=round(hi, 1), ds=[ds_r, ds_d],
                data=q.flatten().tolist())


# --------------------------------------------------------------- server
PAGE = r"""<!doctype html><html><head><meta charset="utf-8">
<title>passiveTuna scope</title>
<style>
:root{--bg:#0e1116;--panel:#151a21;--ink:#d8dde3;--mut:#7e8996;--acc:#e0922f;--data:#56b3c7;--line:#262d36}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 "IBM Plex Mono",Consolas,monospace}
header{display:flex;gap:1.2rem;align-items:baseline;padding:.7rem 1rem;border-bottom:2px solid var(--acc)}
header h1{font:800 1.05rem/1 Archivo,Arial,sans-serif;margin:0;letter-spacing:.02em}
header .b{color:var(--mut);font-size:.75rem}
main{display:grid;grid-template-columns:1fr 22rem;gap:.8rem;padding:.8rem 1rem}
.card{background:var(--panel);border:1px solid var(--line);border-radius:4px;padding:.7rem}
canvas{width:100%;height:auto;background:#000;display:block;border-radius:3px}
.ctl{display:flex;gap:.6rem;align-items:center;margin:.5rem 0;flex-wrap:wrap}
select,button{background:#1c232c;color:var(--ink);border:1px solid var(--line);border-radius:3px;padding:.25rem .5rem;font:inherit}
input[type=range]{flex:1;min-width:12rem}
table{border-collapse:collapse;width:100%;font-size:.78rem}
th{color:var(--mut);text-align:left;font-weight:500;border-bottom:1px solid var(--line);padding:.2rem .3rem}
td{padding:.18rem .3rem;border-bottom:1px solid #1b2128;white-space:nowrap;font-variant-numeric:tabular-nums}
.k{color:var(--mut);font-size:.75rem}
.leg span{display:inline-block;margin-right:.8rem}
.ring{color:#ffd24a}.xh{color:#ff6b6b}.comb{color:#6a7580}.trk{color:var(--data)}
</style></head><body>
<header><h1>passiveTuna · range-Doppler scope</h1><span class="b" id="build"></span><span class="b" id="cap"></span></header>
<main>
<section class="card">
 <div class="ctl">
  <label>run <select id="run"></select></label>
  <button id="play">▶</button>
  <input type="range" id="dwell" min="0" max="0" value="0">
  <span id="tlab" class="k"></span>
 </div>
 <canvas id="c" width="1100" height="620"></canvas>
 <div class="leg k" style="margin-top:.4rem">
  <span class="ring">○ CFAR detection</span><span class="xh">✕ ADS-B aircraft, predicted from its own telemetry</span>
  <span class="comb">┆ pilot comb tooth</span><span class="trk">— kinematic track</span>
  <span id="maplab"></span>
 </div>
 <p class="k" id="note"></p>
</section>
<aside class="card">
 <div class="k" style="margin-bottom:.4rem">aircraft predicted inside the map this dwell</div>
 <table><thead><tr><th>callsign</th><th>R_b km</th><th>f_d Hz</th><th>nearest det</th></tr></thead><tbody id="ac"></tbody></table>
 <div class="k" style="margin:.8rem 0 .4rem">kinematic tracks (whole run)</div>
 <table><thead><tr><th>frames</th><th>km</th><th>obs m/s</th><th>pred m/s</th><th>off</th></tr></thead><tbody id="tr"></tbody></table>
</aside>
</main>
<script>
const $=s=>document.querySelector(s);
let RUN=null, MAP=null, playing=false, timer=null;
const C=$('#c'), X=C.getContext('2d');
const PAD={l:64,r:14,t:14,b:40};
const RMAX=60000, FMAX=1000;
function xr(r){return PAD.l+(r/RMAX)*(C.width-PAD.l-PAD.r)}
function yf(f){return PAD.t+((FMAX-f)/(2*FMAX))*(C.height-PAD.t-PAD.b)}
function heat(v){ // 0..255 -> colour (dark navy -> amber -> white)
 const t=v/255; let r,g,b;
 if(t<.5){const u=t/.5;r=10+u*200;g=14+u*110;b=40+u*10}else{const u=(t-.5)/.5;r=210+u*45;g=124+u*120;b=50+u*190}
 return [r|0,g|0,b|0]}
async function api(u){const r=await fetch(u);return r.json()}
async function loadRuns(){
 const rs=await api('/api/runs'); const sel=$('#run'); sel.innerHTML='';
 for(const r of rs){const o=document.createElement('option');o.value=r.name;o.textContent=r.name+(r.has_maps?'':'  (no maps)');sel.appendChild(o)}
 const q=new URLSearchParams(location.search);
 // default: the URL's ?run=, else the LAST run that has maps, else the last run
 let pick=q.get('run'); if(!pick||!rs.some(r=>r.name===pick)){const withMaps=rs.filter(r=>r.has_maps);pick=(withMaps.length?withMaps[withMaps.length-1]:rs[rs.length-1])?.name}
 if(pick){sel.value=pick;await loadRun(pick);const d=parseInt(q.get('dwell')||'0');if(d>0&&d<=+$('#dwell').max){$('#dwell').value=d;await show(d)}}
}
async function loadRun(name){
 RUN=await api('/api/run?name='+encodeURIComponent(name));
 $('#build').textContent=RUN.build; $('#cap').textContent=RUN.capture;
 $('#dwell').max=RUN.times.length-1; $('#dwell').value=0; $('#note').textContent=RUN.note;
 const tb=$('#tr'); tb.innerHTML='';
 for(const t of RUN.tracks){const tr=document.createElement('tr');tr.innerHTML=`<td>${t.n}</td><td>${(t.range0_m/1e3).toFixed(2)}→${(t.range1_m/1e3).toFixed(2)}</td><td>${t.observed_rate_ms}</td><td>${t.predicted_rate_ms}</td><td>${t.agreement_pct}%</td>`;tb.appendChild(tr)}
 await show(0);
}
async function show(i){
 const t=RUN.times[i]; $('#tlab').textContent=`dwell ${i+1}/${RUN.times.length}  t = ${t.toFixed(3)} s`;
 MAP=RUN.has_maps?await api(`/api/map?name=${encodeURIComponent(RUN.name)}&t=${t}`):null;
 $('#maplab').textContent=MAP?`map ${MAP.lo_db}…${MAP.hi_db} dB over local noise`:'(no saved maps for this run — detections only)';
 draw(t);
}
function draw(t){
 X.fillStyle='#000';X.fillRect(0,0,C.width,C.height);
 const W=C.width-PAD.l-PAD.r, H=C.height-PAD.t-PAD.b;
 if(MAP){ // heatmap: rows=range, cols=doppler; stretch to plot area
  const [nr,nd]=MAP.shape, img=X.createImageData(nd,nr);
  for(let i=0;i<nr*nd;i++){const c=heat(MAP.data[i]);img.data[4*i]=c[0];img.data[4*i+1]=c[1];img.data[4*i+2]=c[2];img.data[4*i+3]=255}
  const off=document.createElement('canvas');off.width=nd;off.height=nr;off.getContext('2d').putImageData(img,0,0);
  X.save();X.translate(PAD.l,PAD.t);
  // rotate so range is X and doppler is Y (+ up): draw transposed
  X.imageSmoothingEnabled=true;
  X.translate(0,H);X.rotate(-Math.PI/2);X.drawImage(off,0,0,nd,nr,0,0,H,W); // now range along the old y -> x
  X.restore();
 }
 // comb teeth
 X.strokeStyle='#3a424c';X.setLineDash([3,5]);
 for(const r of RUN.comb_m){if(r>RMAX)break;X.beginPath();X.moveTo(xr(r),PAD.t);X.lineTo(xr(r),PAD.t+H);X.stroke()}
 X.setLineDash([]);
 // zero-doppler notch
 X.strokeStyle='#2a3139';X.beginPath();X.moveTo(PAD.l,yf(0));X.lineTo(PAD.l+W,yf(0));X.stroke();
 // axes
 X.fillStyle='#9aa4b0';X.font='11px monospace';
 for(let r=0;r<=RMAX;r+=10000){X.fillText((r/1e3)+' km',xr(r)-10,C.height-PAD.b+16)}
 for(let f=-FMAX;f<=FMAX;f+=250){X.fillText((f>0?'+':'')+f+' Hz',4,yf(f)+4)}
 // tracks
 X.strokeStyle='#56b3c7';X.lineWidth=1.5;
 for(const tr of RUN.tracks){X.beginPath();tr.points.forEach((p,k)=>{const x=xr(p.range_m),y=yf(p.doppler_hz);k?X.lineTo(x,y):X.moveTo(x,y)});X.stroke()}
 X.lineWidth=1;
 // detections at this dwell
 const det=RUN.detections.filter(d=>Math.abs(d.t_s-t)<1e-3);
 X.strokeStyle='#ffd24a';
 for(const d of det){X.beginPath();X.arc(xr(d.range_m),yf(d.doppler_hz),5+Math.max(0,(d.snr_db-11))*0.6,0,6.283);X.stroke()}
 // predicted aircraft at this dwell
 const tb=$('#ac');tb.innerHTML='';
 for(const a of RUN.aircraft){
  const row=(RUN.preds[a.icao]||[]).find(p=>Math.abs(p[0]-t)<1e-3); if(!row)continue;
  const x=xr(row[1]),y=yf(row[2]);
  X.strokeStyle='#ff6b6b';X.beginPath();X.moveTo(x-7,y-7);X.lineTo(x+7,y+7);X.moveTo(x-7,y+7);X.lineTo(x+7,y-7);X.stroke();
  X.fillStyle='#ff9b9b';X.fillText(a.callsign,x+9,y-6);
  let best=null;for(const d of det){const dr=Math.abs(d.range_m-row[1]),df=Math.abs(d.doppler_hz-row[2]);const s=dr/43.4+df/3.9;if(!best||s<best.s)best={s,dr,df}}
  const tr=document.createElement('tr');
  const near=best?`${best.dr.toFixed(0)} m / ${best.df.toFixed(0)} Hz`:'—';
  const hit=best&&best.dr<490&&best.df<30;
  tr.innerHTML=`<td>${a.callsign}</td><td>${(row[1]/1e3).toFixed(2)}</td><td>${row[2].toFixed(0)}</td><td style="color:${hit?'#8fe38f':'#7e8996'}">${near}${hit?' ✓':''}</td>`;
  tb.appendChild(tr);
 }
}
$('#dwell').addEventListener('input',e=>show(+e.target.value));
$('#run').addEventListener('change',e=>loadRun(e.target.value));
$('#play').addEventListener('click',()=>{playing=!playing;$('#play').textContent=playing?'❚❚':'▶';
 if(playing){timer=setInterval(()=>{const s=$('#dwell');let v=+s.value+1;if(v>+s.max)v=0;s.value=v;show(v)},280)}else clearInterval(timer)});
loadRuns();
</script></body></html>"""


class H(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype="application/json"):
        self.send_response(code)
        self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body if isinstance(body, bytes) else body.encode())

    def do_GET(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        try:
            if u.path == "/":
                return self._send(200, PAGE, "text/html")
            if u.path == "/api/runs":
                return self._send(200, json.dumps(runs()))
            if u.path == "/api/run":
                run, _ = load_run(q["name"][0])
                return self._send(200, json.dumps(run))
            if u.path == "/api/map":
                m = map_json(q["name"][0], float(q["t"][0]))
                return self._send(200, json.dumps(m))
            return self._send(404, json.dumps({"error": "no such route"}))
        except Exception as e:                                 # noqa: BLE001
            return self._send(500, json.dumps({"error": f"{type(e).__name__}: {e}"}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8650)
    a = ap.parse_args()
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), H)
    print(f"[{BUILD}] http://127.0.0.1:{a.port}/   runs: "
          f"{', '.join(r['name'] for r in runs()) or '(none yet)'}", flush=True)
    srv.serve_forever()


if __name__ == "__main__":
    main()
