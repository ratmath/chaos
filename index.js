// @ratmath/chaos: dynamical systems and chaos for RatMath.
//
// Maps of the line so far: true orbits beside floating-point ones, periodic
// points as exact algebraic numbers, Lyapunov exponents along true orbits,
// bifurcations proven, and symbolic dynamics. `plan` lists what is coming.

export { plan } from "./src/plan.js";
export { ChaosError } from "./src/tree.js";
export { IteratedMap, map, logistic, tent, quadratic, doubling } from "./src/map.js";
export { lyapunov, sensitivity, bifurcationDiagram, bifurcationPoint, superstable, periodDoubling, periodicWindow, itinerary, kneading, sharkovskii } from "./src/dynamics.js";
