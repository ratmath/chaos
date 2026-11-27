# Changelog

## 0.3.0

- Phase 3, maps of the plane: `planeMap`, `henon`, `standardMap` with exact,
  enclosed and float orbits and `compare`; `periodicPoints(n, region)`,
  `fixedPoints` and `cycles` prove each point by Krawczyk's test on Fⁿ − id
  and say whether the list is complete; `exactPeriodicPoints(n)` gives
  algebraic numbers when one coordinate can be eliminated; the Jacobian
  determinant is found symbolically when it is constant. `catMap`, `torusMap`
  and `horseshoe` are exact throughout.
- Phase 3, flows: `flow`, `lorenz`, `rossler`, `vanDerPol`, `pendulum`.
  `solve` encloses the true solution with a Taylor method in outward-rounded
  float intervals that carries the derivative with respect to the starting
  point; `returnMap` encloses a section's return map and its derivative;
  `periodicOrbit` proves a closed orbit exists, with its period, multipliers
  and stability; `equilibria` finds and classifies every equilibrium in a
  box; `bifurcations` brackets the parameters where equilibria appear, merge
  or change stability and names the change.
- `plan` is now empty: every phase is built.
- Docs: five new topic pages and two new explore panels.
- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

## 0.2.0

- Phase 4: iterated function systems of affine maps with rational
  coefficients (`IFS`, `Affine`, `attractors`): exact stages and fixed points,
  rational contraction bounds, the similarity dimension exact when rational
  and certified otherwise (`moran` solves Σ rᵢ^s = 1 for unequal ratios), the
  chaos game with a seeded generator and `boxCounting` (both estimated).
  Cantor sets by allowed digits (`cantor`): exact membership of any fraction,
  stages, measure, dimension. L-systems (`lsystem`, `lsystems`) with symbol
  counts by matrix power and exact turtle coordinates in Q(√3). Elementary
  cellular automata (`elementary`) and Life-like rules (`life`, `patterns`)
  with `classify`.
- Phase 5: `mandelbrot`, `julia` and `newtonBasins` classify every pixel as
  proved inside, proved outside or undecided, with outward-rounded float
  intervals, trapping boxes and the mean-value Newton map;
  `mandelbrotPoint` and `juliaPoint` decide single exact points with the
  complex package's rectangles. Adds `@ratmath/complex` as a vendored
  dependency.
- Docs: six new topic pages and three new explore panels.

## 0.1.0

- Maps of the line (roadmap phase 1): `map`, `logistic`, `tent`, `quadratic`,
  `doubling`. Orbits exact (`orbit`), certified (`enclosedOrbit`, which adds
  working digits until the answer is as accurate as asked) and in doubles
  (`floatOrbit`), with `compare` giving the step at which the doubles part
  from the true orbit. Fixed and periodic points of polynomial maps as exact
  algebraic numbers with multipliers and stability (`fixedPoints`,
  `periodicPoints`, `cycles`); multipliers of exactly 0 and ±1 are recognised
  exactly. Certified periodic points of any formula (`periodicPointsIn`).
- Phase 2: `lyapunov` (the average of log|f′| along the true orbit, enclosed),
  `sensitivity`, `bifurcationPoint` (period doublings and tangent
  bifurcations proven by Krawczyk's test), `periodDoubling` with the ratios
  that approach Feigenbaum's constant, `periodicWindow`, `superstable`,
  `bifurcationDiagram` (estimated), `itinerary`, `kneading`, `sharkovskii`.
- Docs: overview, explore, reference and five topic pages.
- Depends on `@ratmath/analytic` and `@ratmath/algebraic`, vendored as git
  submodules.
- Package skeleton and docs site.
