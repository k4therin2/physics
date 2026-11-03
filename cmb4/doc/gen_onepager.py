import html, json
REPO="https://github.com/k4therin2/physics/tree/cmb-fourth-order"
BLOB="https://github.com/k4therin2/physics/blob/cmb-fourth-order/"
# (id, role, model, tokens_k, tool_uses, minutes, note)
agents=[
("Research","research","Fable 5.1",219,89,10,"Anthropic writeups on physics/math harnesses"),
("Literature","research","Fable 5.1",470,101,17,"33 papers pulled; conventions, kernels, gaps"),
("C1 Einstein eqs","build","Fable 5.1",270,29,31,"sympy engine, Einstein+fluid to 3rd order"),
("C2 Observable","build","Fable 5.1",278,48,48,"shear integral along the ray to 3rd order"),
("C3 Harmonics","build","Fable 5.1",232,34,32,"a_lm / C_l machinery, lensing control"),
("N1 Planar solver","build","Fable 5.1",229,31,33,"exact numerical GR + radiation, ray tracing"),
("W1 Writer","write","Fable 5.1",266,26,22,"settled sections of the manuscript"),
("V1 Referee","verify","Fable 5.1",142,25,26,"hostile referee of the exact formula"),
("V3 Referee","verify","Opus 5.5",231,50,41,"independent check of Einstein engine"),
("V2 Squeezed","verify","Opus 5.5",183,26,71,"adiabatic-mode route to squeezed limits"),
("C5 Linear C_l","build","Opus 5.5",204,42,29,"first-order spectrum, closed-form ISW"),
("C4 Closed forms","build","Opus 5.5",270,106,139,"all 2nd-order fields in closed form"),
("A1 Plane-wave test","assemble","Opus 5.5",241,50,28,"analytic chain vs exact solver (orders 1-3)"),
("A2 K2 + C22","assemble","Opus 5.5",266,65,31,"second-order kernel, C22"),
("A3 K3c + C13","assemble","Opus 5.5",426,117,94,"contracted third-order kernel, C13"),
("V4 Monte Carlo","verify","Opus 5.5",125,30,87,"180k-realisation check of C22 assembly"),
("V5 Zero modes","verify","Opus 5.5",347,112,78,"gauge zero-mode invariance proof"),
("V6 Residues","verify","Opus 5.5",231,36,27,"analytic collinear-pole cancellation"),
("R1 Cold read","referee","Opus 5.5",249,22,13,"manuscript referee #1"),
("W2a Writer","write","Opus 5.5",306,63,17,"map + checks sections"),
("W2b Writer","write","Opus 5.5",281,55,14,"assembly + monochromatic (2nd order)"),
("W3 Revision","write","Opus 5.5",487,105,32,"referee-driven revision, roadmap"),
("W2c Writer","write","Opus 5.5",278,60,14,"third-order + unresolved sections"),
("R2 Cold read","referee","Opus 5.5",334,28,12,"manuscript referee #2"),
("W4 Final edit","write","Opus 5.5",512,141,32,"final editorial pass"),
("Novelty check","research","Opus 5.5",93,16,4,"prior-work search"),
]
roles={"research":("Research","#2a78d6","#3987e5"),"build":("Build","#eb6834","#d95926"),"verify":("Verify","#1baf7a","#199e70"),
       "assemble":("Assemble","#eda100","#c98500"),"write":("Write","#e87ba4","#d55181"),"referee":("Referee","#6250d6","#9085e9")}
tot=sum(a[3] for a in agents); fable=sum(a[3] for a in agents if a[2]=="Fable 5.1"); opus=tot-fable
mins=sum(a[5] for a in agents)
# bar chart svg (horizontal bars, sorted by tokens desc)
srt=sorted(agents,key=lambda a:-a[3])
W=760; left=170; right=70; bh=18; gap=6; top=10
H=top+len(srt)*(bh+gap)+30
maxt=max(a[3] for a in agents); scale=(W-left-right)/maxt
bars=[]
for i,a in enumerate(srt):
    y=top+i*(bh+gap); w=a[3]*scale
    cls="r-"+a[1]
    bars.append(f'<g class="bar" data-name="{html.escape(a[0])}" data-tok="{a[3]}" data-tools="{a[4]}" data-min="{a[5]}" data-model="{a[2]}" data-note="{html.escape(a[6])}"><rect x="{left}" y="{y}" width="{w:.1f}" height="{bh}" rx="0" class="{cls}"/><rect x="{left+w-4:.1f}" y="{y}" width="4" height="{bh}" rx="2" class="{cls}"/><text x="{left-8}" y="{y+bh*0.72:.1f}" class="lbl" text-anchor="end">{html.escape(a[0])}</text><text x="{left+w+6:.1f}" y="{y+bh*0.72:.1f}" class="val">{a[3]}k{" ·O" if a[2].startswith("Opus") else " ·F"}</text></g>')
axis=''.join(f'<line x1="{left+t*scale:.1f}" y1="{top}" x2="{left+t*scale:.1f}" y2="{H-28}" class="grid"/><text x="{left+t*scale:.1f}" y="{H-12}" class="tick" text-anchor="middle">{t}k</text>' for t in (0,100,200,300,400,500))
barsvg=f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" aria-label="Tokens per agent">{axis}{"".join(bars)}</svg>'
legend=''.join(f'<span class="lg"><i class="sw r-{k}"></i>{v[0]}</span>' for k,v in roles.items())
# workflow svg: columns of stages
stages=[("Brief","research",["Problem statement","Research","Literature","Master spec"]),
 ("Build","build",["C1 Einstein eqs","C2 Observable","C3 Harmonics","N1 Planar solver","C4 Closed forms","C5 Linear C_l"]),
 ("Verify","verify",["V1 Referee","V3 Referee","V2 Squeezed"]),
 ("Assemble","assemble",["A1 Plane-wave test","A2 K2 + C22","A3 K3c + C13"]),
 ("Verify","verify",["V4 Monte Carlo","V5 Zero modes","V6 Residues"]),
 ("Write + referee","write",["W1","W2a/b/c","R1","W3","R2","W4"])]
cw=118; ch=26; cg=8; xg=20; sw=len(stages)*(cw+xg); sh=max(len(s[2]) for s in stages)*(ch+cg)+40
els=[]
for i,(name,role,items) in enumerate(stages):
    x=i*(cw+xg)
    els.append(f'<text x="{x+cw/2}" y="14" class="stage" text-anchor="middle">{name}</text>')
    for j,it in enumerate(items):
        y=26+j*(ch+cg); r=role
        if name=="Write + referee" and it.startswith("R"): r="referee"
        els.append(f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" rx="4" class="node r-{r}"/><text x="{x+cw/2}" y="{y+17}" class="nlbl" text-anchor="middle">{html.escape(it)}</text>')
    if i<len(stages)-1:
        els.append(f'<path d="M{x+cw+2},{26+ch/2+30} l{xg-4},0" class="arrow" marker-end="url(#ah)"/>')
wfsvg=f'<svg viewBox="0 0 {sw} {sh}" width="100%" role="img" aria-label="Workflow"><defs><marker id="ah" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 z" class="ahead"/></marker></defs>{"".join(els)}</svg>'
rows=''.join(f'<tr><td>{html.escape(a[0])}</td><td>{roles[a[1]][0]}</td><td>{a[2]}</td><td class="num">{a[3]}k</td><td class="num">{a[4]}</td><td class="num">{a[5]}</td><td>{html.escape(a[6])}</td></tr>' for a in agents)
page=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CMB Fourth Order</title>
<meta name="description" content="One-page report: a multi-agent Claude run on the fourth-order CMB temperature power spectrum problem">
<style>
:root{{--bg:#fcfcfb;--card:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--ink3:#8a8984;--line:#e4e3df;--accent:#2a78d6;
--research:#2a78d6;--build:#eb6834;--verify:#1baf7a;--assemble:#eda100;--write:#e87ba4;--referee:#6250d6}}
@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{--bg:#1a1a19;--card:#232322;--ink:#f0efec;--ink2:#c3c2b7;--ink3:#8f8e87;--line:#36352f;--accent:#3987e5;
--research:#3987e5;--build:#d95926;--verify:#199e70;--assemble:#c98500;--write:#d55181;--referee:#9085e9}}}}
:root[data-theme="dark"]{{--bg:#1a1a19;--card:#232322;--ink:#f0efec;--ink2:#c3c2b7;--ink3:#8f8e87;--line:#36352f;--accent:#3987e5;
--research:#3987e5;--build:#d95926;--verify:#199e70;--assemble:#c98500;--write:#d55181;--referee:#9085e9}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}}
main{{max-width:860px;margin:0 auto;padding:32px 16px 64px}}h1{{font-size:1.7rem;line-height:1.2;margin:0 0 6px}}h2{{font-size:1.15rem;margin:36px 0 10px;border-bottom:1px solid var(--line);padding-bottom:6px}}
p{{margin:0 0 12px}}.sub{{color:var(--ink2);margin-bottom:18px}}a{{color:var(--accent)}}
.tiles{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}}
.tile{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px}}.tile b{{display:block;font-size:1.5rem;line-height:1.1}}.tile span{{color:var(--ink2);font-size:.85rem}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:12px 0}}
ul{{padding-left:20px}}li{{margin:4px 0}}.links a{{display:inline-block;margin:4px 12px 4px 0}}
.r-research{{fill:var(--research)}}.r-build{{fill:var(--build)}}.r-verify{{fill:var(--verify)}}.r-assemble{{fill:var(--assemble)}}.r-write{{fill:var(--write)}}.r-referee{{fill:var(--referee)}}
.node{{stroke:var(--card);stroke-width:2;opacity:.92}}.nlbl{{font-size:10.5px;fill:#fff;font-weight:600}}.stage{{font-size:12px;fill:var(--ink2);font-weight:600}}
.arrow{{stroke:var(--ink3);stroke-width:1.5;fill:none}}.ahead{{fill:var(--ink3)}}
.lbl{{font-size:11px;fill:var(--ink)}}.val{{font-size:10.5px;fill:var(--ink2)}}.tick{{font-size:10px;fill:var(--ink3)}}.grid{{stroke:var(--line);stroke-width:1}}
.bar rect:first-child{{opacity:.88}}.bar:hover rect{{opacity:1}}.lg{{display:inline-flex;align-items:center;gap:6px;margin-right:14px;font-size:.85rem;color:var(--ink2)}}.sw{{width:12px;height:12px;border-radius:3px;display:inline-block}}
.sw.r-research{{background:var(--research)}}.sw.r-build{{background:var(--build)}}.sw.r-verify{{background:var(--verify)}}.sw.r-assemble{{background:var(--assemble)}}.sw.r-write{{background:var(--write)}}.sw.r-referee{{background:var(--referee)}}
#tip{{position:fixed;pointer-events:none;background:var(--card);border:1px solid var(--line);border-radius:6px;padding:8px 10px;font-size:.82rem;color:var(--ink);display:none;max-width:280px;box-shadow:0 4px 14px rgba(0,0,0,.15)}}
table{{width:100%;border-collapse:collapse;font-size:.82rem}}th,td{{text-align:left;padding:5px 6px;border-bottom:1px solid var(--line);vertical-align:top}}.num{{text-align:right;font-variant-numeric:tabular-nums}}
details summary{{cursor:pointer;color:var(--ink2)}}.small{{font-size:.85rem;color:var(--ink2)}}
.pre{{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:10px 12px;font:.82rem/1.45 ui-monospace,Menlo,Consolas,monospace;overflow-x:auto;white-space:pre}}
.kw{{font-weight:600}}
</style></head><body>
<main>
<h1>I threw my Claude at your CMB problem</h1>
<p class="sub">A two-day multi-agent run on “Analytical temperature power spectrum through fourth order” (radiation-only GR universe). This page is the summary; everything behind it is in the repo.</p>
<div class="links card"><span class="kw">Where to look:</span>
<a href="{BLOB}cmb4/doc/SUMMARY.md">2-page summary (SUMMARY.md)</a>
<a href="{BLOB}cmb4/doc/main.pdf">Manuscript PDF (109 pp, roadmap in the first 10)</a>
<a href="{REPO}/cmb4/notes">Agent notes &amp; referee reports</a>
<a href="{REPO}/cmb4/code">Code (sympy engines, solver, checks)</a>
<a href="{REPO}">Branch cmb-fourth-order</a></div>

<h2>Two-sentence takeaway</h2>
<div class="card"><p>For a pure radiation fluid the observed temperature is <em>exactly</em> the local temperature times exp of minus the fluid shear integrated along the ray, which collapses the Sachs–Wolfe / Doppler / lensing / emission-surface bookkeeping into one integral and makes the third-order expansion tractable. The full fourth-order spectrum (C<sub>11</sub>, C<sub>22</sub> and the previously uncomputed 1×3 term) is assembled with all observer and normalisation terms, evaluated for a monochromatic spectrum, and checked end-to-end against an exact numerical GR solution in a different gauge at first, second and third order.</p></div>

<h2>How it ran</h2>
<div class="tiles">
<div class="tile"><b>{len(agents)}</b><span>agent runs (plus 7 that died on a usage limit and were relaunched)</span></div>
<div class="tile"><b>{tot/1000:.1f}M</b><span>agent tokens reported ({fable/1000:.1f}M Fable 5.1, {opus/1000:.1f}M Opus 5.5)</span></div>
<div class="tile"><b>{mins/60:.0f} h</b><span>summed agent wall-clock (many ran in parallel)</span></div>
<div class="tile"><b>~$150–400</b><span>rough API-equivalent cost estimate; it ran on a subscription, so check the usage page</span></div>
</div>
<p>The lead session (Fable 5.1) wrote a master spec and briefs, did no heavy math itself, and kept running notes on disk so a crash cost nothing. Workers ran in six roles. Every claimed result went to a separate referee agent told to assume it was wrong; computer algebra and numerics were the arbiter, never the model’s opinion.</p>
<div class="card">{wfsvg}</div>
<p class="small">Flow left to right; the Verify columns fed corrections back into Build and Assemble. Writers and referees alternated twice.</p>

<h2>Agent size (tokens per run)</h2>
<div class="card">{legend}<div style="margin-top:8px">{barsvg}</div><p class="small">Suffix ·F = Fable 5.1, ·O = Opus 5.5. Hover a bar for tool calls, minutes and the task.</p>
<details><summary>Table view</summary><table><thead><tr><th>Agent</th><th>Role</th><th>Model</th><th class="num">Tokens</th><th class="num">Tool calls</th><th class="num">Minutes</th><th>Task</th></tr></thead><tbody>{rows}</tbody></table></details></div>

<h2>How much of your document got addressed</h2>
<div class="card">
<p><span class="kw">Done:</span> Gaussian variable and “no higher-order primordial modes” defined precisely; initial data satisfying the constraints at every order; emission on T = T<sub>*</sub>, fluid-worldline observer at fixed proper time from a physical initial event, Fermi–Walker tetrad (shown to drop out); normalise-then-remove-dipole in the stated order with all observer and endpoint terms; Θ through third order with every induced vector/tensor term (the 359-term third-order functional is printed in full); every Fourier, convolution and angular normalisation fixed; finite and distinct η<sub>i</sub>, η<sub>*</sub>, η<sub>o</sub>; monochromatic specialisation with the variance held fixed; explicit list of unresolved steps. All analytic; numerics appear only as checks and for one illustrative table.</p>
<p><span class="kw">Not done / partial:</span> third-order time integrals left as explicit one-dimensional integrals (second order is fully closed-form); the 1×3 term evaluated at one parameter point and at η<sub>o</sub>/η<sub>*</sub> = 4 rather than a CMB-like 50; ultraviolet behaviour for power-law spectra only probed numerically (no damping in a perfect fluid, so a log divergence is expected); no test under a general second-order gauge transformation, only the k = 0 family plus the different-gauge numerical comparison.</p></div>

<h2>Is any of it new?</h2>
<div class="card"><ul>
<li><span class="kw">The shear formula:</span> a special case of an exact line-of-sight result in Dunsby 1997 (CQG 14, 3391), who wrote it for general matter and then linearised. The explicit exponential-shear form and its corollary (shear-free radiation flow ⇒ isotropic sky for every comoving observer) were not found stated anywhere; the manuscript credits Dunsby.</li>
<li><span class="kw">The fourth-order spectrum:</span> two literature searches found no complete one-loop CMB temperature spectrum including the 1×3 term (Tomita 2007 explicitly left it out; SONG’s scope excludes it). Appears new, but in a toy radiation-only model.</li>
<li><span class="kw">Reusable tools:</span> a sympy Einstein-equation engine to third order, closed forms for every second-order field in radiation domination, and a high-precision planar-symmetric GR + radiation solver with ray tracing.</li></ul></div>

<h2>Possible next steps (rough API-equivalent guesses)</h2>
<div class="card"><ul>
<li>Evaluate at CMB-like η<sub>o</sub>/η<sub>*</sub> ≈ 50 with a k<sub>0</sub> scan — ~$50–150</li>
<li>Closed forms for the third-order time integrals — ~$100–200</li>
<li>Proper ultraviolet analysis for power-law spectra — ~$50–100</li>
<li>Extend beyond pure radiation (matter + radiation) — a much larger project, likely $500+</li></ul>
<p class="small">Dollar figures are a rough conversion of reported agent tokens and could be off by 2×.</p></div>

<h2>The headline numbers (monochromatic spectrum, k<sub>0</sub>η<sub>*</sub> = 2, k<sub>0</sub>η<sub>o</sub> = 8)</h2>
<div class="pre">  l    C11 (σ²)    ¼C22 (σ⁴)    ⅓Re C13 (σ⁴)    ¼C22 + ⅓Re C13
  2    0.9546      14.2575      −5.3980         +8.8595
  3    0.6251       2.8781     −23.5152        −20.6370
  4    0.1751       1.8659      −9.0639         −7.1980
  6    1.77e-4      0.5221      −1.59e-2        +0.5062
 10    9.68e-6      4.01e-2     −1.46e-3        +3.86e-2</div>
<p class="small">C<sub>ℓ</sub> = ε²C11 + ε⁴[¼C22 + ⅓Re C13]; the ε⁴ term is a correction and can be negative. Fastest way to judge quality: the A1-vs-N1 table in the manuscript’s Checks section (analytic maps vs exact solver, orders 1–3, four configurations).</p>
</main>
<div id="tip"></div>
<script>
const tip=document.getElementById('tip');
document.querySelectorAll('.bar').forEach(b=>{{
 b.addEventListener('mousemove',e=>{{const d=b.dataset;tip.innerHTML='<b>'+d.name+'</b><br>'+d.model+' · '+d.tok+'k tokens · '+d.tools+' tool calls · '+d.min+' min<br>'+d.note;tip.style.display='block';tip.style.left=(e.clientX+14)+'px';tip.style.top=(e.clientY+14)+'px';}});
 b.addEventListener('mouseleave',()=>tip.style.display='none');}});
</script>
</body></html>'''
open('index.html','w').write(page); print(len(page), tot, fable, opus, mins)
