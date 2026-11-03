import html, sys
sys.path.insert(0, '/home/user/physics/cmb4/doc')
OUT = '/home/user/physics/cmb4/doc/onepager_v2.html'
REPO = "https://github.com/k4therin2/physics/tree/cmb-fourth-order"
BLOB = "https://github.com/k4therin2/physics/blob/cmb-fourth-order/"
MAIN = "https://github.com/k4therin2/physics"
E = html.escape

# (id, role, model, tokens_k, tool_uses, minutes, note)  -- copied verbatim from gen_onepager.py
agents = [
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
roles = {"research":"Research","build":"Build","verify":"Verify","assemble":"Assemble","write":"Write","referee":"Referee"}
tot = sum(a[3] for a in agents); fable = sum(a[3] for a in agents if a[2]=="Fable 5.1"); opus = tot-fable
mins = sum(a[5] for a in agents)
nF = sum(1 for a in agents if a[2]=="Fable 5.1"); nO = len(agents)-nF
byname = {a[0]: a for a in agents}

# ---------------------------------------------------------------- workflow SVGs
def node(x, y, w, h, label, role, title=None, fs=12):
    t = f'<title>{E(title)}</title>' if title else ''
    stripe = f'<rect x="{x}" y="{y}" width="4" height="{h}" rx="2" class="f-{role}"/>' if role else ''
    cls = 'wn' if role else 'wn lead-n'
    return (f'<g>{t}<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5" class="{cls}"/>{stripe}'
            f'<text x="{x+12}" y="{y+h/2+fs*0.36:.1f}" class="wt" style="font-size:{fs}px">{E(label)}</text></g>')

def tip(n):
    a = byname.get(n)
    return f'{a[0]} · {a[2]} · {a[3]}k tokens · {a[6]}' if a else None

# Stage definitions: (title, fable nodes, opus nodes); node = (label, role, agent-name-for-title)
STAGES = [
 ("Brief", [("Research","research","Research"),("Literature","research","Literature"),("Master spec","","")], []),
 ("Build in parallel", [("C1 Einstein eqs","build","C1 Einstein eqs"),("C2 Observable","build","C2 Observable"),("C3 Harmonics","build","C3 Harmonics"),("N1 Planar solver","build","N1 Planar solver")],
                       [("C4 Closed forms","build","C4 Closed forms"),("C5 Linear C_l","build","C5 Linear C_l")]),
 ("Referees attack", [("V1 Referee","verify","V1 Referee")], [("V3 Referee","verify","V3 Referee"),("V2 Squeezed","verify","V2 Squeezed")]),
 ("Assemble", [], [("A1 Plane-wave test","assemble","A1 Plane-wave test"),("A2 K2 + C22","assemble","A2 K2 + C22"),("A3 K3c + C13","assemble","A3 K3c + C13")]),
 ("Verify", [], [("V4 Monte Carlo","verify","V4 Monte Carlo"),("V5 Zero modes","verify","V5 Zero modes"),("V6 Residues","verify","V6 Residues")]),
 ("Write ⇄ referee", [("W1 Writer","write","W1 Writer")],
   [("W2a · W2b · W2c","write","W2a Writer"),("R1 Cold read","referee","R1 Cold read"),("W3 Revision","write","W3 Revision"),
    ("R2 Cold read","referee","R2 Cold read"),("W4 Final edit","write","W4 Final edit"),("Novelty check","research","Novelty check")]),
]
DEFS = ('<defs><marker id="ah" viewBox="0 0 8 8" markerWidth="7" markerHeight="7" refX="7" refY="4" orient="auto-start-reverse">'
        '<path d="M0,0.5 L7.5,4 L0,7.5 z" class="ahead"/></marker>'
        '<marker id="ahf" viewBox="0 0 8 8" markerWidth="7" markerHeight="7" refX="7" refY="4" orient="auto-start-reverse">'
        '<path d="M0,0.5 L7.5,4 L0,7.5 z" class="aheadf"/></marker></defs>')

def wf_desktop():
    W = 960; gut = 104; cw = (W-gut)/6; nw = cw-18; nh = 26; ng = 8; pitch = nh+ng
    o = [DEFS]
    # lead band
    o.append(f'<rect x="0" y="0" width="{W}" height="44" rx="8" class="leadband"/>')
    o.append(f'<text x="16" y="19" class="wlead">Lead session · Fable 5.1</text>')
    o.append(f'<text x="16" y="35" class="wsmall">Writes the master spec and every brief, does no heavy maths itself, keeps running notes on disk so a crash costs nothing.</text>')
    ty = 76
    for i, (title, F, O) in enumerate(STAGES):
        x = gut + i*cw
        o.append(f'<line x1="{x+nw/2:.1f}" y1="44" x2="{x+nw/2:.1f}" y2="{ty-16}" class="disp"/>')
        o.append(f'<text x="{x:.1f}" y="{ty}" class="wstage">{i+1}  {E(title)}</text>')
        if i < 5:
            o.append(f'<path d="M{x+nw+2:.1f},{ty-4} h{cw-nw-6:.1f}" class="fwd" marker-end="url(#ahf)"/>')
    fy = 92; frows = max(len(s[1]) for s in STAGES); oy0 = fy + frows*pitch + 30
    orows = max(len(s[2]) for s in STAGES)
    # tier bands
    o.append(f'<rect x="0" y="{fy-8}" width="{W}" height="{frows*pitch+8}" rx="6" class="tierF"/>')
    o.append(f'<text x="12" y="{fy+12}" class="wtier">Fable 5.1</text><text x="12" y="{fy+28}" class="wsmall">first wave</text><text x="12" y="{fy+42}" class="wsmall">{nF} runs · {fable/1000:.1f}M tok</text>')
    o.append(f'<text x="12" y="{oy0+12}" class="wtier">Opus 5.5</text><text x="12" y="{oy0+28}" class="wsmall">after the crash</text><text x="12" y="{oy0+42}" class="wsmall">{nO} runs · {opus/1000:.1f}M tok</text>')
    cy = fy + frows*pitch + 11
    o.append(f'<line x1="0" y1="{cy}" x2="{W}" y2="{cy}" class="crash"/>')
    lab = 'Usage-limit crash: 7 runs died and were relaunched; workers continue on Opus 5.5'
    o.append(f'<rect x="{gut+cw*1.55:.0f}" y="{cy-9}" width="452" height="18" rx="9" class="crashpill"/><text x="{gut+cw*1.55+226:.0f}" y="{cy+4}" class="wcrash" text-anchor="middle">{lab}</text>')
    for i, (title, F, O) in enumerate(STAGES):
        x = gut + i*cw
        for j, (l, r, n) in enumerate(F):
            o.append(node(x, fy+j*pitch, nw, nh, l, r, tip(n) or 'Lead: master spec and briefs'))
        for j, (l, r, n) in enumerate(O):
            o.append(node(x, oy0+j*pitch, nw, nh, l, r, tip(n)))
        if i == 0:
            o.append(f'<path d="M{x+nw/2:.1f},{fy+2*pitch-ng+1} v{ng-2}" class="fwd" marker-end="url(#ahf)"/>')
        if i == 5:
            for j in range(4):
                yy = oy0 + j*pitch + nh
                o.append(f'<path d="M{x+nw/2:.1f},{yy+1} v{ng-2}" class="fwd" marker-end="url(#ahf)"/>')
    # feedback loops below the Opus row
    by = oy0 + orows*pitch + 4
    def arc(i_from, i_to, depth, label):
        x1 = gut + i_from*cw + nw*0.5; x2 = gut + i_to*cw + nw*0.5
        return (f'<path d="M{x1:.1f},{by} C{x1:.1f},{by+depth} {x2:.1f},{by+depth} {x2:.1f},{by}" class="loop" marker-end="url(#ah)"/>'
                f'<text x="{(x1+x2)/2:.1f}" y="{by+depth*0.75+4:.1f}" class="wloop" text-anchor="middle">{label}</text>')
    o.append(arc(2, 1, 34, 'corrections'))
    o.append(arc(4, 3, 34, 'corrections'))
    o.append(arc(4, 1, 86, 'corrections back to build'))
    o.append(arc(5, 4, 34, 'referee findings'))
    H = by + 86
    return f'<svg class="wf wf-d" viewBox="0 0 {W} {H:.0f}" role="img" aria-label="Workflow: a lead session dispatches briefs; builders work in parallel; referees attack; assembly; verification; writers and referees alternate. Fable 5.1 first wave, Opus 5.5 after a usage-limit crash.">{"".join(o)}</svg>'

def wf_mobile():
    W = 336; lm = 20; colw = 148; gap = 12; xF = lm; xO = lm+colw+gap; nh = 25; ng = 6; pitch = nh+ng; fs = 11.5
    o = [DEFS]
    o.append(f'<rect x="0" y="0" width="{W}" height="58" rx="8" class="leadband"/>')
    o.append(f'<text x="12" y="19" class="wlead">Lead session · Fable 5.1</text>')
    o.append(f'<text x="12" y="35" class="wsmall">Writes the spec and every brief, no heavy maths,</text>')
    o.append(f'<text x="12" y="49" class="wsmall">running notes on disk so a crash costs nothing.</text>')
    y = 84
    o.append(f'<text x="{xF}" y="{y}" class="wtier">Fable 5.1</text><text x="{xF}" y="{y+14}" class="wsmall">first wave · {fable/1000:.1f}M tok</text>')
    o.append(f'<text x="{xO}" y="{y}" class="wtier">Opus 5.5</text><text x="{xO}" y="{y+14}" class="wsmall">after crash · {opus/1000:.1f}M tok</text>')
    top = y + 26
    ys = {}
    yy = top
    for i, (title, F, O) in enumerate(STAGES):
        o.append(f'<text x="{lm}" y="{yy+12}" class="wstage">{i+1}  {E(title)}</text>')
        yy += 20
        ys[i] = yy
        for j, (l, r, n) in enumerate(F):
            o.append(node(xF, yy+j*pitch, colw, nh, l, r, tip(n) or 'Lead: master spec and briefs', fs))
        for j, (l, r, n) in enumerate(O):
            o.append(node(xO, yy+j*pitch, colw, nh, l, r, tip(n), fs))
        if i == 5:
            for j in range(4):
                b = yy + j*pitch + nh
                o.append(f'<path d="M{xO+colw/2:.1f},{b+1} v{ng-2}" class="fwd" marker-end="url(#ahf)"/>')
        rows = max(len(F), len(O), 1)
        yy += rows*pitch + 12
    H = yy
    cx = lm + colw + gap/2
    o.insert(1, f'<line x1="{cx}" y1="{top-40}" x2="{cx}" y2="{H-10}" class="crash"/>')
    # feedback arcs at the left margin: stage 3 -> 2, stage 5 -> 4, stage 5 -> 2
    def larc(i_from, i_to, x, lab):
        y1 = ys[i_from] + 10; y2 = ys[i_to] + 10
        return f'<path d="M{lm-2},{y1} C{x},{y1} {x},{y2} {lm-2},{y2}" class="loop" marker-end="url(#ah)"/>'
    o.append(larc(2, 1, 10, ''))
    o.append(larc(4, 3, 10, ''))
    o.append(larc(4, 1, 2, ''))
    return f'<svg class="wf wf-m" viewBox="0 0 {W} {H:.0f}" role="img" aria-label="Workflow, stacked by stage: Fable 5.1 first wave on the left, Opus 5.5 after the usage-limit crash on the right.">{"".join(o)}</svg>'

# ---------------------------------------------------------------- bar chart (HTML, responsive)
SCALE = 550
def bars_html():
    out = []
    for model, label in (("Fable 5.1", "Fable 5.1 · first wave"), ("Opus 5.5", "Opus 5.5 · after the crash")):
        grp = sorted([a for a in agents if a[2] == model], key=lambda a: -a[3])
        sub = sum(a[3] for a in grp)
        out.append(f'<div class="bgrp"><span>{label}</span><span class="bgrp-n">{len(grp)} runs · {sub/1000:.1f}M tokens</span></div>')
        mx, mn = grp[0][0], grp[-1][0]
        for a in grp:
            show = a[0] in (mx, mn)
            val = f'{a[3]}k' if show else ''
            out.append(f'<div class="brow" tabindex="0" data-name="{E(a[0])}" data-role="{roles[a[1]]}" data-tok="{a[3]}" data-tools="{a[4]}" data-min="{a[5]}" data-model="{a[2]}" data-note="{E(a[6])}" aria-label="{E(a[0])}, {roles[a[1]]}, {a[2]}, {a[3]}k tokens">'
                       f'<span class="bl">{E(a[0])}</span><span class="bt"><i class="bar b-{a[1]}" style="width:{a[3]/SCALE*100:.2f}%"></i><em class="bv{" on" if show else ""}" style="left:calc({a[3]/SCALE*100:.2f}% + 6px)">{val}</em></span></div>')
    grid = ''.join(f'<i style="left:{t/SCALE*100:.3f}%"></i>' for t in (0, 100, 200, 300, 400, 500))
    ticks = ''.join(f'<span class="{"odd" if t%200 else ""}" style="left:{t/SCALE*100:.3f}%">{t}k</span>' for t in (0, 100, 200, 300, 400, 500))
    return f'<div class="bars"><div class="bgrid" aria-hidden="true">{grid}</div>{"".join(out)}<div class="baxis" aria-hidden="true"><span class="bl"></span><span class="bt">{ticks}</span></div></div>'

legend = ''.join(f'<span class="lg"><i class="sw b-{k}"></i>{v}</span>' for k, v in roles.items())
rows = ''.join(f'<tr><td>{E(a[0])}</td><td>{roles[a[1]]}</td><td>{a[2]}</td><td class="num">{a[3]}k</td><td class="num">{a[4]}</td><td class="num">{a[5]}</td><td>{E(a[6])}</td></tr>' for a in agents)

# ---------------------------------------------------------------- physics diagram
def physics_svg():
    o = ['<defs><marker id="pa" viewBox="0 0 8 8" markerWidth="7" markerHeight="7" refX="7" refY="4" orient="auto"><path d="M0,0.5 L7.5,4 L0,7.5 z" class="pray-h"/></marker></defs>']
    o.append('<path d="M36,34 C50,66 26,98 38,132 C50,166 26,198 40,226" class="psurf"/>')
    o.append('<text x="6" y="22" class="ptxt">emission surface</text>')
    o.append('<text x="6" y="248" class="psub">where T<tspan dy="3" font-size="9">phys</tspan><tspan dy="-3"> = T</tspan><tspan dy="3" font-size="9">*</tspan></text>')
    o.append('<path d="M40,133 C100,114 160,150 220,128 S290,118 318,124" class="pray" marker-end="url(#pa)"/>')
    for (cx, cy, rx, ry, rot) in ((90, 125, 14, 8, -30), (150, 133, 15, 7.5, 40), (212, 130, 12.5, 10, -65)):
        o.append(f'<circle cx="{cx}" cy="{cy}" r="10.5" class="pref"/>')
        o.append(f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}" transform="rotate({rot} {cx} {cy})" class="pell"/>')
    o.append('<path d="M150,110 L150,96" class="plead"/><text x="170" y="70" class="ptxt" text-anchor="middle">fluid patches the light crosses are sheared:</text>')
    o.append('<text x="170" y="86" class="psub" text-anchor="middle">stretched one way, squeezed the other (σ)</text>')
    o.append('<circle cx="331" cy="124" r="5.5" class="pobs"/>')
    o.append('<text x="354" y="108" class="ptxt" text-anchor="end">observer</text>')
    o.append('<text x="354" y="152" class="psub" text-anchor="end">moves with the fluid,</text>')
    o.append('<text x="354" y="167" class="psub" text-anchor="end">local temperature T<tspan dy="3" font-size="9">o</tspan></text>')
    o.append('<path d="M44,184 q0,10 10,10 h112 q10,0 14,10 q4,-10 14,-10 h112 q10,0 10,-10" class="pbrace"/>')
    o.append('<text x="180" y="224" class="ptxt" text-anchor="middle">S(n) = shear integrated along the ray</text>')
    o.append('<text x="180" y="239" class="psub" text-anchor="middle">emission surface → observer</text>')
    return '<svg class="phys" viewBox="0 0 360 256" role="img" aria-label="Light ray from the emission surface to a comoving observer, crossing sheared fluid patches; the observed temperature is the local temperature times exp of minus the shear integrated along the ray.">' + ''.join(o) + '</svg>'

# ---------------------------------------------------------------- content blocks
OK = '<svg viewBox="0 0 16 16" class="ic" aria-hidden="true"><circle cx="8" cy="8" r="7.5" class="ic-ok"/><path d="M4.5 8.3l2.2 2.2 4.6-4.8" class="ic-mk"/></svg>'
PART = '<svg viewBox="0 0 16 16" class="ic" aria-hidden="true"><circle cx="8" cy="8" r="7.5" class="ic-part"/><path d="M8 0.5a7.5 7.5 0 0 1 0 15z" class="ic-partfill"/></svg>'
NO = '<svg viewBox="0 0 16 16" class="ic" aria-hidden="true"><circle cx="8" cy="8" r="7.5" class="ic-no"/><path d="M5.2 5.2l5.6 5.6M10.8 5.2l-5.6 5.6" class="ic-nomk"/></svg>'
ST = {"done": (OK, "Done"), "part": (PART, "Partial"), "no": (NO, "Not done")}
checklist = [
 ("done", "Gaussian primordial variable and “no higher-order primordial modes” defined precisely", "§5"),
 ("done", "Initial data satisfying the constraints at every order", "§5"),
 ("done", "Emission on T = T<sub>*</sub>; fluid-worldline observer at fixed proper time from a physical initial event; Fermi–Walker tetrad (shown to drop out)", "§4, §8"),
 ("done", "Normalise, then remove the dipole, in the stated order, with all observer and endpoint terms", "§4, §8–9"),
 ("done", "Θ through third order with every induced vector/tensor term; the 359-term third-order functional printed in full", "§8–9, App. C"),
 ("done", "Every Fourier, convolution and angular normalisation fixed", "§2, §10"),
 ("done", "Finite and distinct η<sub>i</sub>, η<sub>*</sub>, η<sub>o</sub>", "§2, §5"),
 ("done", "Monochromatic specialisation with the variance held fixed", "§12"),
 ("done", "Explicit list of unresolved steps", "§14"),
 ("done", "All analytic; numerics appear only as checks and for one illustrative table", "§13"),
 ("part", "Third-order time integrals: left as explicit one-dimensional integrals (second order is fully closed-form)", "§6, App. B"),
 ("part", "The 1×3 term: evaluated at one parameter point, and at η<sub>o</sub>/η<sub>*</sub> = 4 rather than a CMB-like 50", "§12"),
 ("part", "Ultraviolet behaviour for power-law spectra: only probed numerically (no damping in a perfect fluid, so a log divergence is expected)", "§11, §14"),
 ("no", "Test under a general second-order gauge transformation: not done; only the k = 0 family plus the different-gauge numerical comparison", "§13, §14"),
]
def checklist_html():
    r = []
    prev = None
    for s, req, where in checklist:
        ic, lab = ST[s]
        gap = ' first-gap' if prev == 'done' and s != 'done' else ''
        prev = s
        r.append(f'<tr class="st-{s}{gap}"><td class="cst">{ic}<span>{lab}</span></td><td class="creq">{req}<span class="cwhere-m">{where}</span></td><td class="cwhere">{where}</td></tr>')
    return '<table class="check"><thead><tr><th>Status</th><th>Requirement from the brief</th><th>Manuscript</th></tr></thead><tbody>' + ''.join(r) + '</tbody></table>'

verif = [
 ("Exact shear formula", "the result that organises everything",
  ["Hostile referee (V1), told to assume it was wrong", "Kasner spacetime: agrees to 10<sup>−31</sup>", "Symbolically, through third order", "Exact planar GR + radiation solver: residual 1.4×10<sup>−14</sup> at Θ<sub>rms</sub> = 0.08"], None),
 ("Einstein + fluid engine", "the sympy perturbation engine",
  ["Bianchi identity to ε<sup>3</sup>", "Independent re-derivation by a second referee (V3)", "Literature controls reproduced: Kohri–Terada tensor source, Lu et al. vector constraint, Pitrou–Uzan–Bernardeau initial conditions", "Harmonic machinery reproduces the full-sky lensed C<sub>ℓ</sub> (Hu 2000) exactly"], None),
 ("End to end, orders 1–3", "analytic chain vs an exact solver",
  ["Analytic Θ<sub>1</sub>, Θ<sub>2</sub>, Θ<sub>3</sub> maps vs the exact planar solver, single plane wave, four configurations, in a different gauge", "Agreement within the solver’s noise (a<sub>20</sub>[Θ<sub>3</sub>] to 7×10<sup>−8</sup>)"], "A single plane wave probes only axisymmetric parts."),
 ("C<sub>22</sub> assembly", "Monte Carlo (V4)",
  ["180,000 realisations of a 20-mode shell field, K<sub>2</sub> used as a black box", "Agrees with the quadrature within 1.1σ at 0.3–0.7% for ℓ = 1–11 (ℓ = 12: 2.0σ)"], "Tests the assembly, Wick factors and normalisation, not K<sub>2</sub> itself."),
 ("The manuscript", "two cold-read referees",
  ["R1 cold read → W3 referee-driven revision", "R2 cold read → W4 final editorial pass"], None),
]
def verif_html():
    r = []
    for t, sub, items, cav in verif:
        lis = ''.join(f'<li>{OK}<span>{i}</span></li>' for i in items)
        c = f'<p class="vcav">{PART}<span>{cav}</span></p>' if cav else ''
        r.append(f'<li class="vrow"><div class="vhead"><h3>{t}</h3><span>{sub}</span></div><div class="vbody"><ul>{lis}</ul>{c}</div></li>')
    return '<ol class="vchain">' + ''.join(r) + '</ol>'

numrows = [("2","0.9546","14.2575","−5.3980","+8.8595"),("3","0.6251","2.8781","−23.5152","−20.6370"),
           ("4","0.1751","1.8659","−9.0639","−7.1980"),("6","1.77e-4","0.5221","−1.59e-2","+0.5062"),
           ("10","9.68e-6","4.01e-2","−1.46e-3","+3.86e-2")]
def sci(s):
    if 'e' in s:
        m, e = s.split('e'); e = int(e)
        return f'{m}×10<sup>{"−" if e<0 else ""}{abs(e)}</sup>'
    return s
numtable = ('<table class="nums"><thead><tr><th>ℓ</th><th>C<sub>11</sub> <small>(σ²)</small></th><th>¼C<sub>22</sub> <small>(σ⁴)</small></th>'
            '<th>⅓Re C<sub>13</sub> <small>(σ⁴)</small></th><th>¼C<sub>22</sub> + ⅓Re C<sub>13</sub></th></tr></thead><tbody>' +
            ''.join(f'<tr><td>{l}</td><td>{sci(a)}</td><td>{sci(b)}</td><td>{sci(c)}</td><td class="{"neg" if d.startswith("−") else ""}">{sci(d)}</td></tr>' for l, a, b, c, d in numrows) +
            '</tbody></table>')

CSS = r"""
:root{color-scheme:light;
--page:#f8f7f3;--card:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--ink3:#898781;--line:#e1e0d9;--rule:#c3c2b7;--ring:rgba(11,11,11,.10);
--accent:#a64b24;--accent-wash:#f3e7df;--quote:#2b2a27;--tierF:#f1efe8;--leadband:#efede6;
--research:#2a78d6;--build:#eb6834;--verify:#1baf7a;--assemble:#eda100;--write:#e87ba4;--referee:#6250d6;
--good:#0ca30c;--warn:#fab219;--crit:#d03b3b;--neg:#a64b24}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;
--page:#121211;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--ink3:#898781;--line:#2c2c2a;--rule:#383835;--ring:rgba(255,255,255,.10);
--accent:#e3906b;--accent-wash:#2a1f19;--quote:#e9e7e1;--tierF:#20201e;--leadband:#232321;
--research:#3987e5;--build:#d95926;--verify:#199e70;--assemble:#c98500;--write:#d55181;--referee:#9085e9;--neg:#e3906b}}
:root[data-theme="dark"]{color-scheme:dark;
--page:#121211;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--ink3:#898781;--line:#2c2c2a;--rule:#383835;--ring:rgba(255,255,255,.10);
--accent:#e3906b;--accent-wash:#2a1f19;--quote:#e9e7e1;--tierF:#20201e;--leadband:#232321;
--research:#3987e5;--build:#d95926;--verify:#199e70;--assemble:#c98500;--write:#d55181;--referee:#9085e9;--neg:#e3906b}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--page);color:var(--ink);font:17px/1.6 "Source Sans 3",system-ui,-apple-system,"Segoe UI",sans-serif;overflow-x:hidden}
.serif,h1,h2,h3,.quote{font-family:"Newsreader",Georgia,"Times New Roman",serif}
a{color:var(--accent);text-underline-offset:2px;text-decoration-thickness:1px}
a:hover{text-decoration-thickness:2px}
.wrap{max-width:720px;margin:0 auto;padding:0 16px}
.wide{max-width:1000px;margin:0 auto;padding:0 16px}
.inner{max-width:688px;margin:0 auto}
section{padding:56px 0 8px}
.eyebrow{font-size:.78rem;letter-spacing:.09em;text-transform:uppercase;color:var(--ink3);font-weight:600;margin:0 0 10px}
h1{font-weight:500;font-size:clamp(2.1rem,6vw,3.3rem);line-height:1.08;letter-spacing:-.015em;margin:0 0 16px}
h2{font-weight:500;font-size:clamp(1.55rem,3.6vw,2rem);line-height:1.2;letter-spacing:-.01em;margin:0 0 14px}
h3{font-weight:600;font-size:1.12rem;line-height:1.3;margin:0}
p{margin:0 0 14px}
.lede{font-size:1.15rem;color:var(--ink2);margin:0 0 8px}
.lede b{color:var(--ink);font-weight:600}
.secnum{display:block;font:600 .78rem/1 "Source Sans 3",system-ui,sans-serif;letter-spacing:.09em;text-transform:uppercase;color:var(--accent);margin-bottom:10px}
.muted{color:var(--ink2)}.small{font-size:.88rem;color:var(--ink2)}
hr.r{border:0;border-top:1px solid var(--line);margin:0}
/* hero */
header.hero{padding:64px 0 8px}
.quote{position:relative;margin:34px 0 0;padding:6px 0 6px 26px;border-left:3px solid var(--accent);font-size:clamp(1.15rem,2.6vw,1.38rem);line-height:1.5;color:var(--quote);font-weight:400}
.quote p{margin:0 0 12px}.quote p:last-child{margin:0}
.quote cite{display:block;margin-top:12px;font:600 .78rem/1 "Source Sans 3",system-ui,sans-serif;letter-spacing:.09em;text-transform:uppercase;color:var(--ink3);font-style:normal}
.quick{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:34px 0 0}
.quick div{background:var(--card);border:1px solid var(--ring);border-radius:10px;padding:14px 16px;font-size:.95rem;color:var(--ink2)}
.quick b{display:block;color:var(--ink);font-weight:600;margin-bottom:2px}
nav.toc{margin:28px 0 0;font-size:.9rem;color:var(--ink3);display:flex;flex-wrap:wrap;gap:4px 16px}
nav.toc a{color:var(--ink2);text-decoration:none;border-bottom:1px solid var(--line)}
nav.toc a:hover{color:var(--accent);border-color:var(--accent)}
/* cards */
.card{background:var(--card);border:1px solid var(--ring);border-radius:12px;padding:20px}
.fig{margin:22px 0 10px}
.figcap{font-size:.88rem;color:var(--ink2);margin:10px 2px 0}
/* physics */
.physgrid{display:grid;grid-template-columns:1fr;gap:18px;align-items:center}
.phys{width:100%;max-width:470px;height:auto;display:block;margin:0 auto}
.psurf{fill:none;stroke:var(--accent);stroke-width:2.5;stroke-linecap:round}
.pray{fill:none;stroke:var(--ink);stroke-width:1.6}.pray-h{fill:var(--ink)}
.pref{fill:none;stroke:var(--ink3);stroke-width:1;stroke-dasharray:2 3}
.pell{fill:var(--accent-wash);stroke:var(--accent);stroke-width:1.6}
.pobs{fill:var(--ink)}
.plead{stroke:var(--ink3);stroke-width:1}
.pbrace{fill:none;stroke:var(--ink3);stroke-width:1.2}
.ptxt{font:600 12.5px "Source Sans 3",system-ui,sans-serif;fill:var(--ink)}
.psub{font:12px "Source Sans 3",system-ui,sans-serif;fill:var(--ink2)}
.formula{font-family:"Newsreader",Georgia,serif;font-size:clamp(1.15rem,3.4vw,1.45rem);text-align:center;margin:16px 0 4px;color:var(--ink);line-height:1.4}
.formula i{font-style:italic}
.formula small{display:block;font:14px/1.5 "Source Sans 3",system-ui,sans-serif;color:var(--ink2);margin-top:4px}
/* tiles */
.tiles{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:22px 0}
.tile{background:var(--card);border:1px solid var(--ring);border-radius:12px;padding:16px 16px 14px}
.tile .tl{font-size:.85rem;color:var(--ink2);display:block}
.tile .tv{display:block;font:600 2rem/1.15 "Source Sans 3",system-ui,sans-serif;letter-spacing:-.01em;margin:4px 0 6px;color:var(--ink)}
.tile .td{font-size:.82rem;color:var(--ink2);line-height:1.4;display:block}
.seg{display:flex;gap:2px;height:8px;margin:8px 0 6px}
.seg i{display:block;height:8px;border-radius:4px}
.seg .sF{background:var(--ink3)}.seg .sO{background:var(--ink2)}
.segk{display:flex;justify-content:space-between;font-size:.78rem;color:var(--ink2)}
/* workflow */
.wf{width:100%;height:auto;display:block}
.wf-m{display:none;max-width:420px;margin:0 auto}
.leadband{fill:var(--leadband)}
.tierF{fill:var(--tierF)}
.wn{fill:var(--card);stroke:var(--rule);stroke-width:1}
.lead-n{fill:var(--leadband);stroke:var(--ink3);stroke-dasharray:3 2}
.wt{font-family:"Source Sans 3",system-ui,sans-serif;fill:var(--ink)}
.wlead{font:600 13px "Source Sans 3",system-ui,sans-serif;fill:var(--ink)}
.wsmall{font:12px "Source Sans 3",system-ui,sans-serif;fill:var(--ink2)}
.wstage{font:600 12.5px "Source Sans 3",system-ui,sans-serif;fill:var(--ink);letter-spacing:.01em}
.wtier{font:600 13px "Newsreader",Georgia,serif;fill:var(--ink)}
.wcrash{font:600 11.5px "Source Sans 3",system-ui,sans-serif;fill:var(--accent)}
.crash{stroke:var(--accent);stroke-width:1.2;stroke-dasharray:5 4}
.crashpill{fill:var(--card);stroke:var(--accent);stroke-width:1}
.disp{stroke:var(--rule);stroke-width:1}
.fwd{stroke:var(--ink3);stroke-width:1.2;fill:none}.aheadf{fill:var(--ink3)}
.loop{stroke:var(--verify);stroke-width:1.4;fill:none}.ahead{fill:var(--verify)}
.wloop{font:italic 12px "Source Sans 3",system-ui,sans-serif;fill:var(--ink2)}
.f-research{fill:var(--research)}.f-build{fill:var(--build)}.f-verify{fill:var(--verify)}.f-assemble{fill:var(--assemble)}.f-write{fill:var(--write)}.f-referee{fill:var(--referee)}
.wkey{display:flex;flex-wrap:wrap;gap:6px 18px;margin:14px 0 0;font-size:.85rem;color:var(--ink2)}
.wkey span{display:inline-flex;align-items:center;gap:7px}
.k-loop{width:22px;height:0;border-top:1.5px solid var(--verify)}
.k-crash{width:22px;height:0;border-top:1.5px dashed var(--accent)}
/* bars */
.lgs{display:flex;flex-wrap:wrap;gap:6px 16px;margin:0 0 14px}
.lg{display:inline-flex;align-items:center;gap:7px;font-size:.86rem;color:var(--ink2)}
.sw{width:10px;height:10px;border-radius:2px;display:inline-block}
.b-research{background:var(--research)}.b-build{background:var(--build)}.b-verify{background:var(--verify)}.b-assemble{background:var(--assemble)}.b-write{background:var(--write)}.b-referee{background:var(--referee)}
.bars{--lab:150px;--valw:52px;position:relative;font-size:.86rem}
.bgrid{position:absolute;top:0;bottom:22px;left:calc(var(--lab) + 12px);right:var(--valw);pointer-events:none}
.bgrid i{position:absolute;top:0;bottom:0;width:0;border-left:1px solid var(--line)}
.bgrid i:first-child{border-left-color:var(--rule)}
.bgrp{position:relative;display:flex;justify-content:space-between;align-items:baseline;gap:10px;margin:14px 0 4px;padding:0 0 4px;font-weight:600;font-size:.88rem;color:var(--ink);background:linear-gradient(var(--card),var(--card)) left/ calc(var(--lab) + 12px) 100% no-repeat}
.bgrp:first-of-type{margin-top:0}
.bgrp-n{font-weight:400;color:var(--ink2);background:var(--card);padding-left:8px}
.brow{display:grid;grid-template-columns:var(--lab) 1fr;column-gap:12px;align-items:center;height:24px;border-radius:4px;outline:none;padding-right:var(--valw)}
.brow:hover,.brow:focus-visible{background:color-mix(in srgb,var(--ink) 5%,transparent)}
.brow:focus-visible{box-shadow:0 0 0 2px var(--accent)}
.bl{text-align:right;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:var(--ink)}
.bt{position:relative;height:12px}
.bar{position:absolute;left:0;top:0;height:12px;border-radius:0 4px 4px 0;display:block}
.bv{position:absolute;top:-3px;font-style:normal;font-size:.8rem;color:var(--ink2);font-variant-numeric:tabular-nums;white-space:nowrap}
.baxis{display:grid;grid-template-columns:var(--lab) 1fr;column-gap:12px;height:22px;margin-top:4px;padding-right:var(--valw)}
.baxis .bt{height:22px}
.baxis .bt span{position:absolute;top:2px;transform:translateX(-50%);font-size:.76rem;color:var(--ink3);font-variant-numeric:tabular-nums}
#tip{position:fixed;z-index:10;pointer-events:none;background:var(--card);border:1px solid var(--ring);border-radius:8px;padding:9px 11px;font-size:.84rem;line-height:1.45;color:var(--ink);display:none;max-width:290px;box-shadow:0 6px 20px rgba(0,0,0,.14)}
#tip b{font-weight:600}#tip .tm{color:var(--ink2)}
details{margin-top:16px}
details summary{cursor:pointer;color:var(--ink2);font-size:.9rem}
.tscroll{overflow-x:auto;-webkit-overflow-scrolling:touch;margin-top:10px}
table{border-collapse:collapse;width:100%}
.data{font-size:.82rem;min-width:620px}
.data th,.data td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}
.data th{font-weight:600;color:var(--ink2)}
.num{text-align:right!important;font-variant-numeric:tabular-nums}
/* checklist */
.check{font-size:.95rem}
.check th{font:600 .74rem/1 "Source Sans 3",system-ui,sans-serif;letter-spacing:.08em;text-transform:uppercase;color:var(--ink3);text-align:left;padding:0 10px 10px;border-bottom:1px solid var(--rule)}
.check td{padding:11px 10px;border-bottom:1px solid var(--line);vertical-align:top}
.check tr.first-gap td{border-top:1.5px solid var(--rule)}
.cst{white-space:nowrap;width:112px;font-weight:600;font-size:.86rem}
.cst span{vertical-align:middle}
.cst .ic{margin-right:7px}
.st-done .cst{color:var(--ink)}
.creq{color:var(--ink)}
.st-part .creq,.st-no .creq{color:var(--ink)}
.cwhere{white-space:nowrap;color:var(--ink2);font-size:.86rem;width:110px;font-variant-numeric:tabular-nums}
.cwhere-m{display:none}
.ic{width:16px;height:16px;vertical-align:-3px;flex:none}
.ic-ok{fill:var(--good)}.ic-mk{fill:none;stroke:#fff;stroke-width:1.8;stroke-linecap:round;stroke-linejoin:round}
.ic-part{fill:none;stroke:#c98500;stroke-width:1.4}.ic-partfill{fill:#c98500}
.ic-no{fill:var(--crit)}.ic-nomk{stroke:#fff;stroke-width:1.8;stroke-linecap:round}
.tally{display:flex;flex-wrap:wrap;gap:6px 18px;margin:0 0 16px;font-size:.92rem;color:var(--ink2)}
.tally span{display:inline-flex;align-items:center;gap:7px}
/* novelty */
.cards3{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin-top:20px}
.ncard{background:var(--card);border:1px solid var(--ring);border-radius:12px;padding:20px;display:flex;flex-direction:column}
.ncard .tag{align-self:flex-start;font:600 .72rem/1 "Source Sans 3",system-ui,sans-serif;letter-spacing:.07em;text-transform:uppercase;padding:5px 8px;border-radius:999px;margin-bottom:12px;border:1px solid var(--rule);color:var(--ink2)}
.ncard .tag.t-new{border-color:var(--accent);color:var(--accent)}
.ncard h3{margin-bottom:8px}
.ncard p{font-size:.95rem;color:var(--ink2);margin:0}
.ncard p b{color:var(--ink);font-weight:600}
/* verification chain */
.vchain{list-style:none;margin:20px 0 0;padding:0;position:relative}
.vrow{position:relative;display:grid;grid-template-columns:220px 1fr;gap:20px;padding:18px 0 18px 30px;border-bottom:1px solid var(--line)}
.vrow:last-child{border-bottom:0}
.vrow::before{content:"";position:absolute;left:7px;top:0;bottom:0;border-left:1.5px solid var(--rule)}
.vrow:first-child::before{top:26px}.vrow:last-child::before{bottom:calc(100% - 26px)}
.vrow::after{content:"";position:absolute;left:2px;top:21px;width:9px;height:9px;border-radius:50%;background:var(--card);border:2px solid var(--verify)}
.vhead span{display:block;font-size:.86rem;color:var(--ink2);margin-top:2px}
.vbody ul{list-style:none;margin:0;padding:0}
.vbody li{display:flex;gap:9px;align-items:flex-start;margin:0 0 6px;font-size:.95rem}
.vbody li .ic,.vcav .ic{margin-top:4px}
.vcav{display:flex;gap:9px;align-items:flex-start;font-size:.9rem;color:var(--ink2);margin:8px 0 0}
/* numbers */
.nums{font-size:.95rem;font-variant-numeric:tabular-nums;min-width:520px}
.nums th{font-weight:600;font-size:.86rem;color:var(--ink2);text-align:right;padding:0 12px 10px;border-bottom:1px solid var(--rule);white-space:nowrap}
.nums td{text-align:right;padding:9px 12px;border-bottom:1px solid var(--line);white-space:nowrap}
.nums th:first-child,.nums td:first-child{text-align:left;padding-left:2px}
.nums td:last-child{font-weight:600}
.nums td.neg{color:var(--neg)}
.nums small{font-weight:400;color:var(--ink3)}
.callout{background:var(--accent-wash);border-radius:12px;padding:16px 18px;margin:18px 0 0;font-size:.98rem}
.callout b{font-weight:600}
/* next steps */
.steps{list-style:none;margin:16px 0 0;padding:0}
.steps li{display:flex;justify-content:space-between;align-items:baseline;gap:16px;padding:12px 0;border-bottom:1px solid var(--line)}
.steps li span:last-child{white-space:nowrap;font-weight:600;font-variant-numeric:tabular-nums}
/* links */
.links{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px;margin-top:18px}
.links a{display:block;background:var(--card);border:1px solid var(--ring);border-radius:12px;padding:14px 16px;text-decoration:none;color:var(--ink)}
.links a:hover{border-color:var(--accent)}
.links a b{display:block;font-weight:600;color:var(--accent)}
.links a span{font-size:.88rem;color:var(--ink2)}
footer{margin:64px 0 0;padding:24px 0 48px;border-top:1px solid var(--line);font-size:.85rem;color:var(--ink3)}
@media (max-width:860px){
 .tiles{grid-template-columns:repeat(2,1fr)}
 .cards3{grid-template-columns:1fr}
 .vrow{grid-template-columns:1fr;gap:8px}
}
@media (max-width:640px){
 body{font-size:16px}
 section{padding:44px 0 4px}
 header.hero{padding:40px 0 4px}
 .quote{padding-left:18px}
 .card{padding:16px}
 .wf-d{display:none}.wf-m{display:block}
 .wfcard{padding:12px 6px}
 .wfcard .wkey{padding:0 8px}
 .bars{--lab:118px;--valw:40px;font-size:.8rem}
 .cwhere{display:none}
 .cwhere-m{display:block;font-size:.82rem;color:var(--ink2);margin-top:3px}
 .cst{width:auto;padding-right:4px}
 .cst span{display:none}
 .check th:first-child{font-size:0;padding:0}
 .check th:last-child{display:none}
 .seg + .segk{flex-direction:column;gap:0}
 .baxis .bt span.odd{display:none}
 .nums{min-width:0;font-size:.8rem}
 .nums th{white-space:normal;padding:0 5px 8px;font-size:.74rem;vertical-align:bottom}
 .nums td{padding:8px 5px}
 .numcard{padding:12px 10px}
 .tile .tv{font-size:1.7rem}
 .steps li{flex-direction:column;gap:2px}
}
"""

JS = r"""
const tip=document.getElementById('tip');
function show(b,x,y){const d=b.dataset;tip.innerHTML='<b>'+d.name+'</b> <span class="tm">· '+d.role+'</span><br>'+d.model+' · '+d.tok+'k tokens · '+d.tools+' tool calls · '+d.min+' min<br><span class="tm">'+d.note+'</span>';tip.style.display='block';
 const w=tip.offsetWidth,h=tip.offsetHeight;let L=x+14,T=y+14;if(L+w>innerWidth-8)L=Math.max(8,x-w-14);if(T+h>innerHeight-8)T=Math.max(8,y-h-14);tip.style.left=L+'px';tip.style.top=T+'px';}
document.querySelectorAll('.brow').forEach(b=>{
 b.addEventListener('mousemove',e=>show(b,e.clientX,e.clientY));
 b.addEventListener('mouseleave',()=>tip.style.display='none');
 b.addEventListener('focus',()=>{const r=b.getBoundingClientRect();show(b,r.left+r.width*0.4,r.bottom-8);});
 b.addEventListener('blur',()=>tip.style.display='none');});
"""

page = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CMB Fourth Order</title>
<meta name="description" content="One-page report: a two-day multi-agent Claude run on the fourth-order CMB temperature power spectrum in a radiation-only universe — what was done, how, what is verified, what is new, and what is not done.">
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,500;0,6..72,600;1,6..72,400&family=Source+Sans+3:ital,wght@0,400;0,600;1,400&display=swap" rel="stylesheet">
<style>{CSS}</style></head><body>

<header class="hero"><div class="wrap">
<p class="eyebrow">I threw my Claude at your CMB problem</p>
<h1>The CMB temperature spectrum through fourth order</h1>
<p class="lede">A <b>two-day multi-agent Claude run</b> on a fourth-order CMB problem: “Analytical temperature power spectrum through fourth order”, in a radiation-only GR universe. This page is the summary; everything behind it is in the repo.</p>
<blockquote class="quote">
<p>For a pure radiation fluid the observed temperature is <em>exactly</em> the local temperature times exp of minus the fluid shear integrated along the ray, which collapses the Sachs–Wolfe / Doppler / lensing / emission-surface bookkeeping into one integral and makes the third-order expansion tractable.</p>
<p>The full fourth-order spectrum (C<sub>11</sub>, C<sub>22</sub> and the previously uncomputed 1×3 term) is assembled with all observer and normalisation terms, evaluated for a monochromatic spectrum, and checked end-to-end against an exact numerical GR solution in a different gauge at first, second and third order.</p>
<cite>The two-sentence takeaway</cite>
</blockquote>
<div class="quick">
<div><b>If you have five minutes</b>Open the A1-vs-N1 table in the manuscript’s Checks section (§13): analytic maps vs an exact solver, orders 1–3, four configurations.</div>
<div><b>What’s not done</b>Third-order time integrals not closed; 1×3 term at one parameter point; UV only probed numerically; no general gauge test. <a href="#coverage">Details</a></div>
</div>
<nav class="toc" aria-label="On this page"><a href="#problem">The problem</a><a href="#how">How it ran</a><a href="#size">Agent size</a><a href="#coverage">Coverage</a><a href="#new">Is it new?</a><a href="#verify">Verification</a><a href="#numbers">Numbers</a><a href="#next">Next steps</a><a href="#links">Links</a></nav>
</div></header>

<section id="problem"><div class="wrap">
<span class="secnum">01 · The physics</span>
<h2>What the problem is</h2>
<p>The brief asks for the CMB temperature power spectrum C<sub>ℓ</sub> (how strongly the sky’s temperature varies at each angular scale ℓ) computed analytically through fourth order in the primordial perturbation, in a toy universe that contains only radiation and obeys full general relativity. At fourth order that means not just the textbook first-order term C<sub>11</sub>, but the product of two second-order maps (C<sub>22</sub>) and a first-order map times a third-order one (C<sub>13</sub>, the “1×3” term), with every observer and normalisation term kept.</p>
<p>The usual route tracks Sachs–Wolfe, Doppler, lensing and emission-surface effects separately, and that bookkeeping is what makes third order hard. The run found that for a pure radiation fluid there is an exact shortcut: the temperature seen in direction <i>n</i> is the local fluid temperature at the observer times the exponential of minus the fluid’s shear, integrated along the light ray back to the emission surface. Sachs–Wolfe, ISW, Doppler, lensing, time delay and the location of the emission surface all live inside that one integral, so the observer’s own temperature drops out of the normalised map and the third-order expansion becomes tractable.</p>
<div class="card fig">
{physics_svg()}
<p class="formula"><i>T</i><sub>obs</sub>(<i>n</i>) = <i>T</i><sub>o</sub> · exp[ −<i>S</i>(<i>n</i>) ], &nbsp; <i>S</i>(<i>n</i>) = ∫<sub>e</sub><sup>o</sup> <i>σ</i><sub>μν</sub> <i>e</i><sup>μ</sup><i>e</i><sup>ν</sup> d<i>l</i>
<small>σ = fluid shear, e = photon direction in the fluid frame, dl = proper length in that frame. Normalised map: 1 + Θ = e<sup>−S</sup>/⟨e<sup>−S</sup>⟩<sub>sky</sub>.</small></p>
</div>
<p class="small">The emission surface enters only as the endpoint of the integral; the observer’s frame enters only as a rotation of the sky, which cannot change C<sub>ℓ</sub>.</p>
</div></section>

<section id="how"><div class="wrap">
<span class="secnum">02 · The workflow</span>
<h2>How it ran</h2>
<p>The lead session (Fable 5.1) wrote a master spec and briefs, did no heavy maths itself, and kept running notes on disk so a crash cost nothing. Workers ran in six roles. Every claimed result went to a separate referee agent told to assume it was wrong; computer algebra and numerics were the arbiter, never the model’s opinion.</p>
</div>
<div class="wide"><div class="tiles">
<div class="tile"><span class="tl">Agent runs</span><span class="tv">{len(agents)}</span><span class="td">plus 7 that died on a usage limit and were relaunched</span></div>
<div class="tile"><span class="tl">Agent tokens reported</span><span class="tv">{tot/1000:.1f}M</span><div class="seg" role="img" aria-label="{fable/1000:.1f}M Fable 5.1, {opus/1000:.1f}M Opus 5.5"><i class="sF" style="flex:{fable}"></i><i class="sO" style="flex:{opus}"></i></div><span class="segk"><span>{fable/1000:.1f}M Fable 5.1</span><span>{opus/1000:.1f}M Opus 5.5</span></span></div>
<div class="tile"><span class="tl">Summed agent wall-clock</span><span class="tv">{mins/60:.0f} h</span><span class="td">many ran in parallel</span></div>
<div class="tile"><span class="tl">Rough API-equivalent cost</span><span class="tv">~$150–400</span><span class="td">an estimate: it ran on a subscription, so check the usage page</span></div>
</div>
<div class="card fig wfcard">
{wf_desktop()}
{wf_mobile()}
<div class="wkey"><span class="lg"><i class="sw b-research"></i>Research</span><span class="lg"><i class="sw b-build"></i>Build</span><span class="lg"><i class="sw b-verify"></i>Verify</span><span class="lg"><i class="sw b-assemble"></i>Assemble</span><span class="lg"><i class="sw b-write"></i>Write</span><span class="lg"><i class="sw b-referee"></i>Referee</span><span><i class="k-loop"></i>feedback: corrections back upstream</span><span><i class="k-crash"></i>usage-limit crash</span></div>
</div>
<div class="inner"><p class="figcap">Read left to right (top to bottom on a phone). The verify stages fed corrections back into Build and Assemble; writers and referees alternated twice. Hover a box for its token count and task.</p></div>
</div></section>

<section id="size"><div class="wrap">
<span class="secnum">03 · Agent size</span>
<h2>Tokens per agent run</h2>
<p class="muted">Each bar is one agent run, coloured by role and grouped by model. The two largest runs were the final editing and revision passes.</p>
<div class="card fig">
<div class="lgs">{legend}</div>
{bars_html()}
<p class="figcap">Hover or tab to a bar for tool calls, minutes and the task. Labelled values: largest and smallest run in each group.</p>
<details><summary>Table view (all {len(agents)} runs)</summary><div class="tscroll"><table class="data"><thead><tr><th>Agent</th><th>Role</th><th>Model</th><th class="num">Tokens</th><th class="num">Tool calls</th><th class="num">Minutes</th><th>Task</th></tr></thead><tbody>{rows}</tbody></table></div></details>
</div>
</div></section>

<section id="coverage"><div class="wrap">
<span class="secnum">04 · Against the brief</span>
<h2>Coverage of the original brief</h2>
<div class="tally"><span>{OK}{sum(1 for c in checklist if c[0]=="done")} done</span><span>{PART}{sum(1 for c in checklist if c[0]=="part")} partial</span><span>{NO}{sum(1 for c in checklist if c[0]=="no")} not done</span></div>
<div class="card">{checklist_html()}</div>
<p class="figcap">Section numbers refer to the manuscript PDF (§ = section, App. = appendix).</p>
</div></section>

<section id="new"><div class="wide">
<div class="inner"><span class="secnum">05 · Novelty</span><h2>Is any of it new?</h2>
<p class="muted">Honest answer: partly. Two literature searches were run; here is what they found.</p></div>
<div class="cards3">
<div class="ncard"><span class="tag">Known in general form</span><h3>The shear formula</h3><p>A special case of an exact line-of-sight result in <b>Dunsby 1997</b> (CQG 14, 3391), who wrote it for general matter and then linearised. The explicit exponential-shear form and its corollary (shear-free radiation flow ⇒ isotropic sky for every comoving observer) were not found stated anywhere. The manuscript credits Dunsby.</p></div>
<div class="ncard"><span class="tag t-new">Appears new · toy model</span><h3>The fourth-order spectrum</h3><p>Two literature searches found no complete one-loop CMB temperature spectrum including the 1×3 term (Tomita 2007 explicitly left it out; SONG’s scope excludes it). <b>Appears new, but in a toy radiation-only model.</b></p></div>
<div class="ncard"><span class="tag">Reusable</span><h3>The tools</h3><p>A sympy Einstein-equation engine to third order; closed forms for every second-order field in radiation domination; a high-precision planar-symmetric GR + radiation solver with ray tracing.</p></div>
</div></div></section>

<section id="verify"><div class="wrap">
<span class="secnum">06 · Checks</span>
<h2>Verification at a glance</h2>
<p class="muted">Each link in the chain was checked by something other than the agent that produced it.</p>
{verif_html()}
</div></section>

<section id="numbers"><div class="wrap">
<span class="secnum">07 · The result</span>
<h2>The headline numbers</h2>
<p class="muted">Monochromatic spectrum, k<sub>0</sub>η<sub>*</sub> = 2, k<sub>0</sub>η<sub>o</sub> = 8, variance σ² held fixed.</p>
<p class="formula"><i>C</i><sub>ℓ</sub> = ε²<i>C</i><sub>11</sub> + ε⁴ [ ¼<i>C</i><sub>22</sub> + ⅓Re <i>C</i><sub>13</sub> ]</p>
<div class="card fig numcard"><div class="tscroll" style="margin:0">{numtable}</div></div>
<p class="small">The ε⁴ term is a correction and can be negative (shown in the accent colour).</p>
<div class="callout"><b>Fastest way to judge quality:</b> the A1-vs-N1 table in the manuscript’s Checks section (analytic maps vs exact solver, orders 1–3, four configurations).</div>
</div></section>

<section id="next"><div class="wrap">
<span class="secnum">08 · What next</span>
<h2>Possible next steps</h2>
<p class="muted">Rough API-equivalent guesses.</p>
<ul class="steps">
<li><span>Evaluate at CMB-like η<sub>o</sub>/η<sub>*</sub> ≈ 50 with a k<sub>0</sub> scan</span><span>~$50–150</span></li>
<li><span>Closed forms for the third-order time integrals</span><span>~$100–200</span></li>
<li><span>Proper ultraviolet analysis for power-law spectra</span><span>~$50–100</span></li>
<li><span>Extend beyond pure radiation (matter + radiation): a much larger project</span><span>likely $500+</span></li>
</ul>
<p class="small" style="margin-top:12px">Dollar figures are a rough conversion of reported agent tokens and could be off by 2×.</p>
</div></section>

<section id="links"><div class="wrap">
<span class="secnum">09 · Where to look</span>
<h2>Everything behind this page</h2>
<div class="links">
<a href="{BLOB}cmb4/doc/main.pdf"><b>Manuscript PDF</b><span>109 pp; the roadmap is in the first 10</span></a>
<a href="{BLOB}cmb4/doc/SUMMARY.md"><b>2-page summary</b><span>SUMMARY.md: results, bookkeeping, checks, limitations</span></a>
<a href="{REPO}/cmb4/notes"><b>Agent notes &amp; referee reports</b><span>running notes from the run and the referee reports</span></a>
<a href="{REPO}/cmb4/code"><b>Code</b><span>sympy engines, planar solver, checks</span></a>
<a href="{REPO}"><b>Branch cmb-fourth-order</b><span>the working branch</span></a>
<a href="{MAIN}"><b>Repository (main)</b><span>github.com/k4therin2/physics</span></a>
</div>
<footer>Generated by a Claude Code session. Numbers are agent tokens as reported by the agents; dollar costs are estimates.</footer>
</div></section>

<div id="tip" role="tooltip"></div>
<script>{JS}</script>
</body></html>'''
open(OUT, 'w').write(page)
print(len(page.encode()), tot, fable, opus, mins)
