# Changelog

## Unreleased

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
