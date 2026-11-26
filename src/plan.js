// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Phases 1, 2, 4 and 5 are built; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
  {
    "phase": 3,
    "contents": "Maps of the plane; flows, sections, periodic orbits, bifurcations",
    "estimate": "1,300"
  },
  {
    "phase": 6,
    "contents": "Pictures",
    "estimate": "1,000 in graphics"
  }
].map(Object.freeze));
