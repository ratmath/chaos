// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Nothing here is built yet; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
  {
    "phase": 1,
    "contents": "Maps of the line: orbits, periodic points, stability",
    "estimate": "900"
  },
  {
    "phase": 2,
    "contents": "Bifurcation diagrams, Lyapunov exponents, symbolic dynamics",
    "estimate": "1,000"
  },
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
