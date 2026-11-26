// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Phases 1 and 2 (maps of the line) are built; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
  {
    "phase": 3,
    "contents": "Maps of the plane; flows, sections, periodic orbits, bifurcations",
    "estimate": "1,300"
  },
  {
    "phase": 4,
    "contents": "Iterated function systems, Cantor sets, dimensions, L-systems, cellular automata",
    "estimate": "900"
  },
  {
    "phase": 5,
    "contents": "Julia and Mandelbrot sets, Newton basins",
    "estimate": "600"
  },
  {
    "phase": 6,
    "contents": "Pictures",
    "estimate": "1,000 in graphics"
  }
].map(Object.freeze));
