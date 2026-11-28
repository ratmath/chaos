// @ratmath/chaos: dynamical systems and chaos for RatMath.
//
// Maps of the line: true orbits beside floating-point ones, periodic points as
// exact algebraic numbers, Lyapunov exponents along true orbits, bifurcations
// proven, and symbolic dynamics. Maps of the plane and flows: periodic points
// and periodic orbits proven to exist, equilibria and their bifurcations
// located, the cat map and the horseshoe exactly. Fractals: iterated function systems, Cantor
// sets, dimensions, L-systems and cellular automata; Julia and Mandelbrot sets
// and Newton basins with every pixel classed inside, outside or undecided.
// Every planned phase is built, so `plan` is empty.

export { plan } from "./src/plan.js";
export { ChaosError } from "./src/tree.js";
export { IteratedMap, map, logistic, tent, quadratic, doubling } from "./src/map.js";
export { PlaneMap, planeMap, henon, standardMap, TorusMap, torusMap, catMap, Horseshoe, horseshoe } from "./src/plane.js";
export { Flow, flow, lorenz, rossler, vanDerPol, pendulum } from "./src/flow.js";
export { lyapunov, sensitivity, bifurcationDiagram, bifurcationPoint, superstable, periodDoubling, periodicWindow, itinerary, kneading, sharkovskii } from "./src/dynamics.js";
export { lyapunovExponent, invariantDensity } from "./src/density.js";
export { Affine, IFS, attractors, moran, DigitCantor, cantor, boxCounting, seeded } from "./src/fractals.js";
export { LSystem, lsystem, lsystems, Q3 } from "./src/lsystem.js";
export { Elementary, elementary, showRows, Life, life, patterns } from "./src/automata.js";
export { PixelGrid, mandelbrot, julia, newtonBasins, mandelbrotPoint, juliaPoint } from "./src/escape.js";
