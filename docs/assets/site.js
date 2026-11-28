// Shared page furniture for the chaos docs: the contents column, the theme
// switch and the code boxes. Add each topic page to PAGES as it is written.

import * as A from "./chaos.js?v=b4132d101e";
import { mountRunboxes } from "./runbox.js?v=b4132d101e";

export { A };

export const PAGES = [
  { href: "index.html", title: "Overview" },
  { href: "reference.html", title: "Reference" },
  { href: "explore.html", title: "Explore" },
  { group: "Maps of the line" },
  { href: "orbits.html", title: "Orbits, true and floating", blurb: "Exact and enclosed orbits, and the step where doubles part from them" },
  { href: "periodic-points.html", title: "Periodic points and stability", blurb: "Algebraic periodic points, multipliers, cycles" },
  { href: "lyapunov.html", title: "Lyapunov exponents and sensitivity", blurb: "Certified averages along true orbits; the limit from the invariant density; separation measured" },
  { href: "bifurcations.html", title: "Period doubling and windows", blurb: "Bifurcation parameters proven; Feigenbaum's ratios" },
  { href: "symbolic-dynamics.html", title: "Symbolic dynamics", blurb: "Itineraries, kneading sequences, Sharkovskii's order" },
  { group: "Maps of the plane" },
  { href: "plane-maps.html", title: "Maps of the plane", blurb: "Hénon and standard maps: true orbits, periodic points proven, exact where algebra allows" },
  { href: "cat-map-and-horseshoe.html", title: "The cat map and the horseshoe", blurb: "Chaos done exactly: periods on the torus, symbolic dynamics of Smale's horseshoe" },
  { group: "Flows" },
  { href: "flows.html", title: "Flows, enclosed", blurb: "Solutions of differential equations inside proven boxes; Lorenz, Rössler, Van der Pol" },
  { href: "periodic-orbits.html", title: "Sections and periodic orbits", blurb: "Return maps with derivatives; closed orbits proven to exist, with period and multipliers" },
  { href: "equilibria.html", title: "Equilibria and their bifurcations", blurb: "Every rest point in a box, classified; saddle-node, pitchfork and Hopf parameters bracketed" },
  { group: "Fractals and automata" },
  { href: "iterated-function-systems.html", title: "Iterated function systems", blurb: "Exact stages, contraction bounds, similarity dimensions, the chaos game" },
  { href: "cantor-sets.html", title: "Cantor sets and dimension", blurb: "Exact membership of fractions; stages, measure, dimension" },
  { href: "l-systems.html", title: "L-systems", blurb: "Symbol counts by matrix power; turtle paths exact in Q(√3)" },
  { href: "cellular-automata.html", title: "Cellular automata", blurb: "Elementary rules and Life; still lifes, oscillators, spaceships" },
  { group: "Escape-time pictures" },
  { href: "mandelbrot-and-julia.html", title: "Mandelbrot and Julia sets", blurb: "Pixels proved inside or outside, the rest honestly undecided" },
  { href: "newton-basins.html", title: "Newton basins", blurb: "Each root's basin proved by interval Newton traps" },
];

export const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat(Infinity)) if (kid !== null && kid !== undefined && kid !== false) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};

// ---- chrome -------------------------------------------------------------------

function theme() {
  let saved = null;
  try {
    saved = localStorage.getItem("ratmath-theme");
  } catch {
    /* storage may be blocked */
  }
  const apply = (t) => {
    if (t === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", t);
    try {
      localStorage.setItem("ratmath-theme", t);
    } catch {
      /* ignore */
    }
    for (const b of document.querySelectorAll(".theme button")) b.classList.toggle("on", b.dataset.t === t);
  };
  if (saved) apply(saved);
  return h("div", { class: "theme", role: "group", "aria-label": "Theme" },
    ["auto", "light", "dark"].map((t) => h("button", { "data-t": t, class: (saved ?? "auto") === t ? "on" : "", onclick: () => apply(t) }, t)));
}

function chrome() {
  const here = location.pathname.split("/").pop() || "index.html";
  const aside = document.querySelector("aside");
  if (!aside) return;
  const nav = h("nav", { class: "toc", "aria-label": "Contents" });
  const sel = h("select", { id: "toc-select", "aria-label": "Page", onchange: (e) => (location.href = e.target.value) });
  let group = null;
  for (const p of PAGES) {
    if (p.group) {
      nav.append(h("div", { class: "group" }, p.group));
      group = h("optgroup", { label: p.group });
      sel.append(group);
      continue;
    }
    nav.append(h("a", { href: p.href, class: p.href === here ? "on" : null, "aria-current": p.href === here ? "page" : null }, p.title));
    (group ?? sel).append(h("option", { value: p.href, selected: p.href === here }, p.title));
  }
  aside.append(
    h("header", { class: "brand" }, h("a", { href: "index.html" }, h("div", { class: "mark" }, "RatMath ", h("span", {}, "chaos"))), h("small", {}, "Dynamical systems and chaos")),
    sel, nav, theme(),
  );
}

/** Topic cards for the overview page (pages that have a blurb). */
export function topicCards(el) {
  el.append(...PAGES.filter((p) => p.blurb).map((p) => h("a", { href: p.href }, h("strong", {}, p.title), h("span", {}, p.blurb))));
}

chrome();
mountRunboxes(document, { scope: A });
