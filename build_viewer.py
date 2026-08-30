"""Build a standalone HTML viewer for the three released corpora.

Reads `runs/corpora.json` (see `export_corpus.py`) and embeds it in a single
self-contained page: every formula rendered as the sequence of nodes the
meta-law counting convention actually charges it for, plus the rank ladder each
corpus's `beta` is fitted to.

The corpus is CC0, but it is not vendored into this repo -- run
`export_corpus.py` first, then this. Output lands in `runs/` and is gitignored.

Usage:  .venv/bin/python build_viewer.py
"""
from __future__ import annotations

import json
from pathlib import Path

import export_corpus as X

# Micro-examples answering "how is 2 / pi / sqrt(2) counted?" -- built with the
# same parser as the corpus so they can never drift from it.
EXAMPLES = [
    ("2*pi", "a plain integer and a named constant are two different leaves"),
    ("sqrt(2)", "no single node for the number sqrt(2) -- an operator over a leaf"),
    ("sqrt(2/3)", "and a fraction is a division over two leaves"),
    ("x**2", "an exponent of 2, 3 or 4 folds into the operator; the 2 vanishes"),
    ("x**5", "any other exponent stays a leaf of its own"),
    ("-x", "unary minus is an operator, not part of the number"),
]

TEMPLATE = r"""<title>Formula Node Ledger</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
:root{
  --ground:#F2F3EE; --surface:#FBFBF8; --sunk:#E8EAE2;
  --ink:#191C19; --ink-2:#565D56; --ink-3:#848B84;
  --rule:#D5D9CF; --rule-2:#C2C8BC;
  --accent:#125E66; --accent-soft:#DCEAEA;
  --variable:#2E6B4F; --number:#8A6112; --constant:#6B4A9C;
  --unary:#A8452A; --binary:#2A5D8F; --unknown:#7C837C;
  --flag:#9A5B12; --flag-bg:#F6E9D6;
  --shadow:0 1px 2px rgba(25,28,25,.06),0 4px 14px rgba(25,28,25,.05);
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --ground:#131612; --surface:#1B1F1A; --sunk:#101310;
  --ink:#E9EBE4; --ink-2:#A7AEA6; --ink-3:#7D857C;
  --rule:#2C312A; --rule-2:#3A4038;
  --accent:#5CB9C2; --accent-soft:#17302F;
  --variable:#6FBF92; --number:#D8A64A; --constant:#B292E2;
  --unary:#E5866A; --binary:#7BAEE0; --unknown:#8D948C;
  --flag:#E0A85C; --flag-bg:#2C2415;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 4px 14px rgba(0,0,0,.25);
}}
:root[data-theme="dark"]{
  --ground:#131612; --surface:#1B1F1A; --sunk:#101310;
  --ink:#E9EBE4; --ink-2:#A7AEA6; --ink-3:#7D857C;
  --rule:#2C312A; --rule-2:#3A4038;
  --accent:#5CB9C2; --accent-soft:#17302F;
  --variable:#6FBF92; --number:#D8A64A; --constant:#B292E2;
  --unary:#E5866A; --binary:#7BAEE0; --unknown:#8D948C;
  --flag:#E0A85C; --flag-bg:#2C2415;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 4px 14px rgba(0,0,0,.25);
}
*{box-sizing:border-box}
body{background:var(--ground);color:var(--ink);
  font-family:"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased}
.wrap{max-width:1120px;margin:0 auto;padding:0 24px 96px}
h1,h2,h3{font-family:Spectral,Georgia,"Times New Roman",serif;font-weight:600;
  text-wrap:balance;margin:0}
a{color:var(--accent)}
code,.mono{font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace}
.tnum{font-variant-numeric:tabular-nums}

/* ---------- masthead ---------- */
header.mast{border-bottom:1px solid var(--rule);margin-bottom:34px;
  padding:52px 0 30px}
.eyebrow{font-size:11.5px;letter-spacing:.13em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600;margin-bottom:14px}
h1{font-size:clamp(30px,4.6vw,44px);line-height:1.12;letter-spacing:-.015em}
.standfirst{font-family:Spectral,Georgia,serif;font-size:18.5px;line-height:1.6;
  color:var(--ink-2);max-width:64ch;margin:16px 0 0}
.standfirst em{color:var(--ink)}

/* ---------- stat strip ---------- */
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(168px,1fr));
  gap:1px;background:var(--rule);border:1px solid var(--rule);
  border-radius:3px;overflow:hidden;margin-top:30px}
.stat{background:var(--surface);padding:15px 17px}
.stat .k{font-size:11px;letter-spacing:.1em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600}
.stat .v{font-size:25px;font-weight:500;margin-top:5px;letter-spacing:-.02em}
.stat .n{font-size:12.5px;color:var(--ink-2);margin-top:2px}

/* ---------- sections ---------- */
section{margin-top:52px}
.sec-head{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap;
  border-bottom:1px solid var(--rule);padding-bottom:9px;margin-bottom:20px}
.sec-head h2{font-size:22px;letter-spacing:-.01em}
.sec-head .note{font-size:13px;color:var(--ink-3)}
.prose{max-width:66ch;color:var(--ink-2);font-size:14.5px}
.prose strong{color:var(--ink);font-weight:600}
.prose p{margin:0 0 11px}

/* ---------- node chips ---------- */
.tape{display:flex;flex-wrap:wrap;gap:4px;align-items:center}
.chip{display:inline-flex;align-items:center;gap:5px;font-family:"IBM Plex Mono",monospace;
  font-size:12px;font-weight:500;padding:2.5px 7px;border-radius:3px;
  border:1px solid;line-height:1.45;white-space:nowrap}
.chip .idx{font-size:9.5px;opacity:.62;font-variant-numeric:tabular-nums}
.c-variable{color:var(--variable);border-color:var(--variable);
  background:color-mix(in srgb,var(--variable) 12%,transparent)}
.c-number{color:var(--number);border-color:var(--number);
  background:color-mix(in srgb,var(--number) 12%,transparent)}
.c-constant{color:var(--constant);border-color:var(--constant);
  background:color-mix(in srgb,var(--constant) 12%,transparent)}
.c-unary{color:var(--unary);border-color:var(--unary);
  background:color-mix(in srgb,var(--unary) 12%,transparent)}
.c-binary{color:var(--binary);border-color:var(--binary);
  background:color-mix(in srgb,var(--binary) 12%,transparent)}
.c-unknown{color:var(--unknown);border-color:var(--unknown);border-style:dashed;
  background:color-mix(in srgb,var(--unknown) 9%,transparent)}

/* ---------- worked examples ---------- */
.egs{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:1px;
  background:var(--rule);border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.eg{background:var(--surface);padding:15px 17px;display:flex;flex-direction:column;gap:9px}
.eg .src{font-family:"IBM Plex Mono",monospace;font-size:15px;font-weight:500}
.eg .why{font-size:12.5px;color:var(--ink-2);line-height:1.5}
.eg .cnt{font-size:11px;letter-spacing:.08em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600}

/* ---------- legend ---------- */
.legend{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:18px}
.legend .chip{font-family:"IBM Plex Sans",sans-serif;font-size:12px}

/* ---------- controls ---------- */
.controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;
  position:sticky;top:0;z-index:5;padding:12px 0;background:var(--ground);
  border-bottom:1px solid var(--rule);margin-bottom:18px}
.tabs{display:flex;border:1px solid var(--rule-2);border-radius:3px;overflow:hidden}
.tabs button{font:inherit;font-size:13px;font-weight:500;padding:6px 14px;border:0;
  background:var(--surface);color:var(--ink-2);cursor:pointer;
  border-right:1px solid var(--rule-2)}
.tabs button:last-child{border-right:0}
.tabs button[aria-pressed="true"]{background:var(--accent);color:var(--ground)}
.tabs button:focus-visible,input:focus-visible,select:focus-visible,
.card > summary:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
input[type="search"],select{font:inherit;font-size:13px;padding:6px 10px;
  border:1px solid var(--rule-2);border-radius:3px;background:var(--surface);
  color:var(--ink)}
input[type="search"]{min-width:190px}
.chk{display:inline-flex;align-items:center;gap:6px;font-size:13px;color:var(--ink-2)}
.count{margin-left:auto;font-size:12.5px;color:var(--ink-3);font-variant-numeric:tabular-nums}

/* ---------- formula cards ---------- */
.cards{display:flex;flex-direction:column;gap:9px}
.card{background:var(--surface);border:1px solid var(--rule);border-radius:3px;
  box-shadow:var(--shadow);scroll-margin-top:68px}
.card > summary{cursor:pointer;padding:13px 16px;display:flex;gap:13px;
  align-items:baseline;flex-wrap:wrap;list-style:none}
.card > summary::-webkit-details-marker{display:none}
.card > summary::before{content:"▸";color:var(--ink-3);font-size:11px;
  transition:transform .15s ease;display:inline-block}
@media (prefers-reduced-motion:reduce){.card > summary::before{transition:none}}
.card[open] > summary::before{transform:rotate(90deg)}
.card .id{font-family:"IBM Plex Mono",monospace;font-size:13px;font-weight:600}
.card .corp{font-size:10.5px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600}
.card .kbadge{margin-left:auto;font-size:12px;color:var(--ink-2);
  font-variant-numeric:tabular-nums;white-space:nowrap}
.card .kbadge b{font-weight:600;color:var(--ink);font-size:14px}
.expr{font-family:"IBM Plex Mono",monospace;font-size:13.5px;
  padding:0 16px 13px;color:var(--ink);overflow-x:auto;white-space:pre}
.body{border-top:1px solid var(--rule);padding:15px 16px;display:flex;
  flex-direction:column;gap:16px}
.sub{font-size:11px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600;margin-bottom:8px}
.flagrow{display:flex;flex-wrap:wrap;gap:7px}
.flag{font-size:11.5px;padding:2px 8px;border-radius:3px;background:var(--flag-bg);
  color:var(--flag);border:1px solid color-mix(in srgb,var(--flag) 35%,transparent)}
.tally{display:flex;flex-wrap:wrap;gap:18px;font-size:13px;color:var(--ink-2)}
.tally b{color:var(--ink);font-variant-numeric:tabular-nums}
.tree{font-family:"IBM Plex Mono",monospace;font-size:12.5px;overflow-x:auto}
.tree ul{list-style:none;margin:0;padding-left:17px;border-left:1px solid var(--rule-2)}
.tree > ul{padding-left:0;border-left:0}
.tree li{margin:2px 0}

/* ---------- rank ladder ---------- */
.ladders{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:18px}
.ladder{background:var(--surface);border:1px solid var(--rule);border-radius:3px;
  padding:16px 18px;box-shadow:var(--shadow)}
.ladder h3{font-size:16px;margin-bottom:3px}
.ladder .beta{font-size:12.5px;color:var(--ink-2);margin-bottom:14px}
.ladder .beta b{color:var(--accent);font-weight:600;font-size:14px}
.rung{display:grid;grid-template-columns:18px 74px 1fr 42px;gap:9px;align-items:center;
  margin-bottom:5px;font-size:12px}
.rung .r{color:var(--ink-3);font-variant-numeric:tabular-nums;text-align:right}
.rung .n{font-variant-numeric:tabular-nums;color:var(--ink-2);text-align:right}
.rung .key{font-family:"IBM Plex Mono",monospace;font-size:11px;color:var(--ink-2);
  overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.track{background:var(--sunk);border-radius:2px;overflow:hidden}
.bar{height:15px;border-radius:2px;min-width:2px}
.b-variable{background:var(--variable)} .b-number{background:var(--number)}
.b-constant{background:var(--constant)} .b-unary{background:var(--unary)}
.b-binary{background:var(--binary)}
footer{margin-top:64px;padding-top:18px;border-top:1px solid var(--rule);
  font-size:12.5px;color:var(--ink-3);max-width:70ch}
</style>

<div class="wrap">
<header class="mast">
  <div class="eyebrow">Constantin, Bartlett, Desmond &amp; Ferreira 2024 &middot; arXiv:2408.11065</div>
  <h1>Formula Node Ledger</h1>
  <p class="standfirst">Every formula in the three released corpora, broken into the nodes
  the meta-law counting convention actually charges it for. <em>Complexity is just the number
  of chips.</em> The colour of each chip is the rank bucket it falls into &mdash; and those
  buckets, not the individual operators, are what <span class="mono">beta</span> is fitted to.</p>
  <div class="stats" id="stats"></div>
</header>

<section>
  <div class="sec-head"><h2>How a number becomes nodes</h2>
    <span class="note">the parser, on six deliberately awkward cases</span></div>
  <div class="prose">
    <p>There is no node for &ldquo;the number &radic;2&rdquo; and no node for &ldquo;2&pi;&rdquo;.
    Every numeric idea is built out of the same leaves and operators as everything else, so a
    constant that looks atomic to a physicist can cost four nodes. The one exception is a small
    integer exponent, which is folded into the operator and disappears.</p>
  </div>
  <div class="egs" id="egs"></div>
</section>

<section>
  <div class="sec-head"><h2>The corpora</h2>
    <span class="note">211 formulas &middot; click any row to open its tree</span></div>
  <div class="legend" id="legend"></div>
  <div class="controls">
    <div class="tabs" id="tabs"></div>
    <input type="search" id="q" placeholder="Search name or formula" aria-label="Search formulas">
    <select id="sort" aria-label="Sort order">
      <option value="src">Source order</option>
      <option value="k-desc">Complexity, high to low</option>
      <option value="k-asc">Complexity, low to high</option>
    </select>
    <label class="chk"><input type="checkbox" id="flagged"> Only flagged</label>
    <span class="count" id="count"></span>
  </div>
  <div class="cards" id="cards"></div>
</section>

<section>
  <div class="sec-head"><h2>What beta is actually fitted to</h2>
    <span class="note">operator instances, ranked by frequency</span></div>
  <div class="prose">
    <p>This is the step that surprises people. Every variable in a corpus &mdash; <span
    class="mono">x</span>, <span class="mono">rho</span>, <span class="mono">theta</span>,
    all of them &mdash; collapses into a <strong>single</strong> rank entry called
    <span class="mono">variable</span>. Same for every plain number, and for every named
    constant. So the top of the ladder is dominated by three bucket categories, and
    <span class="mono">beta</span> measures how fast this ladder falls away.</p>
  </div>
  <div class="ladders" id="ladders"></div>
</section>

<footer id="foot"></footer>
</div>

<script type="application/json" id="payload">__DATA__</script>
<script>
const DATA = JSON.parse(document.getElementById("payload").textContent);
const CORPORA = Object.keys(DATA.corpora);
const LABEL = {inflation:"Encyclopaedia Inflationaris", feynman:"Feynman Lectures",
               wiki:"Wikipedia named equations"};
const el = (t,c,x) => { const n=document.createElement(t); if(c)n.className=c;
  if(x!==undefined)n.textContent=x; return n; };

/* Which rank bucket a parsed node falls into. Nullary symbols are only
   resolvable when the row's own tallies pin them down -- see `unknown`. */
function cat(node, roles){
  if(node.arity==="binary") return "binary";
  if(node.arity==="unary")  return "unary";
  if(node.key==="number")   return "number";
  return roles[node.text] || "unknown";
}
function flatten(node, roles, out){
  out = out || [];
  out.push({label: node.text || node.key, cat: cat(node, roles)});
  (node.children||[]).forEach(c => flatten(c, roles, out));
  return out;
}
function tape(nodes, numbered){
  const d = el("div","tape");
  nodes.forEach((n,i) => {
    const c = el("span","chip c-"+n.cat);
    if(numbered) c.appendChild(el("span","idx", i+1));
    c.appendChild(el("span",null,n.label));
    c.title = n.cat==="unknown"
      ? "the CSV gives only a count, not which symbol is the constant"
      : n.cat;
    d.appendChild(c);
  });
  return d;
}

/* ---- stats ---- */
const statBox = document.getElementById("stats");
CORPORA.forEach(c => {
  const d = DATA.corpora[c], ks = d.formulas.map(f=>f.complexity.declared).sort((a,b)=>a-b);
  const s = el("div","stat");
  s.appendChild(el("div","k", LABEL[c]));
  s.appendChild(el("div","v tnum", "beta " + d.beta.toFixed(3)));
  s.appendChild(el("div","n", `${d.formulas.length} formulas · median complexity ${ks[Math.floor(ks.length/2)]}`));
  statBox.appendChild(s);
});

/* ---- worked examples ---- */
const egBox = document.getElementById("egs");
DATA.examples.forEach(e => {
  const nodes = flatten(e.tree, e.roles);
  const d = el("div","eg");
  d.appendChild(el("div","src", e.expr));
  d.appendChild(tape(nodes,false));
  d.appendChild(el("div","cnt", nodes.length + (nodes.length===1?" node":" nodes")));
  d.appendChild(el("div","why", e.why));
  egBox.appendChild(d);
});

/* ---- legend ---- */
const LEG = [["variable","variable — every one, one rank"],["number","plain number"],
  ["constant","named constant (pi, c, …)"],["unary","unary operator"],
  ["binary","binary operator"],["unknown","symbol the CSV leaves ambiguous"]];
const legBox = document.getElementById("legend");
LEG.forEach(([c,t]) => legBox.appendChild(el("span","chip c-"+c, t)));

/* ---- controls ---- */
let state = {corpus:"all", q:"", sort:"src", flagged:false};
const tabBox = document.getElementById("tabs");
[["all","All"], ...CORPORA.map(c=>[c, LABEL[c].split(" ")[0]])].forEach(([v,t])=>{
  const b = el("button",null,t);
  b.setAttribute("aria-pressed", v===state.corpus);
  b.onclick = () => { state.corpus=v;
    [...tabBox.children].forEach((x,i)=>x.setAttribute("aria-pressed", x===b));
    render(); };
  tabBox.appendChild(b);
});
document.getElementById("q").oninput = e => { state.q=e.target.value.toLowerCase(); render(); };
document.getElementById("sort").onchange = e => { state.sort=e.target.value; render(); };
document.getElementById("flagged").onchange = e => { state.flagged=e.target.checked; render(); };

/* ---- formula cards ---- */
const ALL = CORPORA.flatMap(c => DATA.corpora[c].formulas);
const FLAGTEXT = {
  bracket_repaired:"brackets repaired to parse",
  complexity_mismatch:"declared complexity disagrees with the recount",
  number_count_mismatch:"‘numbers’ column disagrees with the literals present"};

function treeEl(node, roles){
  const ul = el("ul");
  const li = el("li");
  const c = cat(node, roles);
  const s = el("span","chip c-"+c, node.text || node.key);
  li.appendChild(s);
  (node.children||[]).forEach(ch => li.appendChild(treeEl(ch, roles)));
  ul.appendChild(li);
  return ul;
}

function card(f){
  const d = el("details","card");
  const sum = el("summary");
  sum.appendChild(el("span","id", f.name));
  sum.appendChild(el("span","corp", f.corpus));
  if(f.flags.length) sum.appendChild(el("span","flag", f.flags.length===1?"1 flag":f.flags.length+" flags"));
  const kb = el("span","kbadge");
  kb.innerHTML = "complexity <b>" + f.complexity.declared + "</b>";
  sum.appendChild(kb);
  d.appendChild(sum);
  d.appendChild(el("div","expr", f.expr));

  const body = el("div","body");
  const nodes = flatten(f.tree, f.symbol_roles);

  const t1 = el("div");
  t1.appendChild(el("div","sub", `The ${nodes.length} counted nodes`));
  t1.appendChild(tape(nodes,true));
  body.appendChild(t1);

  const t2 = el("div");
  t2.appendChild(el("div","sub","Nullary tallies, as the authors filed them"));
  const tal = el("div","tally");
  tal.innerHTML = `<span>variables <b>${f.nullary.variables}</b></span>
    <span>numbers <b>${f.nullary.numbers}</b></span>
    <span>constants <b>${f.nullary.constants}</b></span>
    <span>recount <b>${f.complexity.recomputed}</b> vs declared <b>${f.complexity.declared}</b></span>`;
  t2.appendChild(tal);
  body.appendChild(t2);

  if(f.flags.length){
    const t3 = el("div");
    t3.appendChild(el("div","sub","Flags"));
    const fr = el("div","flagrow");
    f.flags.forEach(x => fr.appendChild(el("span","flag", FLAGTEXT[x]||x)));
    if(f.expr_original) fr.appendChild(el("span","flag","source read: "+f.expr_original));
    t3.appendChild(fr);
    body.appendChild(t3);
  }

  const t4 = el("div");
  t4.appendChild(el("div","sub","Expression tree"));
  const tw = el("div","tree");
  tw.appendChild(treeEl(f.tree, f.symbol_roles));
  t4.appendChild(tw);
  body.appendChild(t4);

  d.appendChild(body);
  return d;
}

function render(){
  const box = document.getElementById("cards");
  box.textContent = "";
  let rows = ALL.filter(f =>
    (state.corpus==="all" || f.corpus===state.corpus) &&
    (!state.flagged || f.flags.length) &&
    (!state.q || f.name.toLowerCase().includes(state.q) || f.expr.toLowerCase().includes(state.q)));
  if(state.sort==="k-desc") rows = [...rows].sort((a,b)=>b.complexity.declared-a.complexity.declared);
  if(state.sort==="k-asc")  rows = [...rows].sort((a,b)=>a.complexity.declared-b.complexity.declared);
  const frag = document.createDocumentFragment();
  rows.forEach(f => frag.appendChild(card(f)));
  box.appendChild(frag);
  document.getElementById("count").textContent =
    rows.length + " of " + ALL.length + " formulas";
}
render();

/* ---- rank ladders ---- */
const ladBox = document.getElementById("ladders");
CORPORA.forEach(c => {
  const d = DATA.corpora[c];
  const box = el("div","ladder");
  box.appendChild(el("h3", null, LABEL[c]));
  const b = el("div","beta");
  b.innerHTML = `fitted <b>beta = ${d.beta.toFixed(3)}</b> &plusmn; ${d.beta_sigma.toFixed(3)}
    over ${d.ranks.length} ranks`;
  box.appendChild(b);
  const max = d.ranks[0].count;
  d.ranks.slice(0,12).forEach((r,i) => {
    const row = el("div","rung");
    row.appendChild(el("div","r", i+1));
    const kind = ["variable","number","constant"].includes(r.key) ? r.key
               : (["Add","Subtract","Mult","Div","Pow"].includes(r.key) ? "binary" : "unary");
    const key = el("div","key", r.key); key.title = r.key;
    row.appendChild(key);
    const track = el("div","track");
    const bar = el("div","bar b-"+kind);
    bar.style.width = Math.max(2, 100*r.count/max) + "%";
    track.appendChild(bar);
    row.appendChild(track);
    row.appendChild(el("div","n", r.count));
    box.appendChild(row);
  });
  ladBox.appendChild(box);
});

document.getElementById("foot").textContent =
  "Data: " + DATA.source + ". Counting convention: " + DATA.convention +
  " Nullary symbols are shown as ambiguous wherever the released columns give only a" +
  " count and not which symbol carries which role.";
</script>
"""


def main() -> None:
    src = Path("runs/corpora.json")
    payload = json.loads(src.read_text())

    payload["examples"] = [
        {"expr": e, "why": why, "tree": X.to_tree(e),
         # pi is the only symbol in these hand-picked cases whose role is
         # unambiguous; everything else is a stand-in variable.
         "roles": {n: ("constant" if n == "pi" else "variable")
                   for n in ["pi", "x"]}}
        for e, why in EXAMPLES
    ]

    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    out = Path("runs/formula-ledger.html")
    out.write_text(TEMPLATE.replace("__DATA__", blob))
    print(f"wrote {out}  ({out.stat().st_size/1024:.0f} kB)")


if __name__ == "__main__":
    main()
