# Changelog

## Unreleased

- Docs: the browser bundle is rebuilt against Core 0.12 (the checked-in one
  still rounded 65 to 60 at the tens and returned a bare BigInt from
  `floor()`), is built without identifier minification (Bun renames classes
  even with `keepNames`, which broke `new Rational(q.floor())` and
  `q.add(q.floor())` in the browser), and the docs build now runs every code
  box against the emitted bundle after a mixed Integer/Rational check. Found
  by the external QA review.

- Requires Core 0.12: the `@ratmath/core` peer range is now `>=0.12.0` (Core 0.12 returns Integers from `floor`/`ceil`/`trunc`/`round` and rounds half up by default).
- Fix: `invariantDensity` and `lyapunovExponent` recognise the logistic map at
  r = 4, the tent map at slope 2 and x² − 2 by their formula, not by the map's
  display name: `map("x", { name: "logistic", params: { r: 4 } })` (the
  identity, exponent 0) was given log 2. The same formula under another name
  or variable is now recognised.
- Neutral periodic points (multiplier exactly ±1) carry `nonlinear` and
  `nonlinearReason`: the Schwarzian derivative decides at −1, the second and
  third derivatives at +1. The logistic fixed point 2/3 at r = 3 is
  "attracting".
- `periodDoubling` refuses with a `budget` message naming the levels it found
  when an engine's BigInt cap is reached, instead of a raw RangeError; the
  bifurcations page runs 5 levels (6 takes seconds, and fails on Bun 1.4.0).
- Docs: the bifurcations page states the nondegeneracy conditions behind the
  ±1 scenarios (+1 also occurs in transcritical and pitchfork bifurcations).
- `lyapunovExponent(f, { distribution, interval })` (Course Test roadmap, item
  B8): the exponent of almost every orbit as a certified enclosure of
  ∫ log|f′| dμ against the invariant measure, given by its distribution
  function. The unbounded density at the ends and the logarithmic singularity
  at a critical point are both bounded rigorously. For `logistic(4)`,
  `tent(2)` and `quadratic(-2)` the measure is built in (`invariantDensity`)
  and the result encloses log 2; a distribution you supply is tested for
  invariance from certified preimages and refused on a mismatch.

- `enclosedOrbit` (and so `lyapunov`, `compare` and the rest that follow one
  exact orbit) took a long numerator for a large value: a bounded orbit whose
  terms carry thousands of bits, such as `lyapunov(logistic(4), "1/5",
  { steps: 2000 })`, was refused as "leaves every bounded range". The test
  now compares the numerator's bits with the denominator's. Found by the
  Course Test (ratmath/shelf#7).

## 0.3.0

- Types: the result of `Flow.bifurcations` is now
  `EquilibriumBifurcation`. It was a second `Bifurcation` interface that
  TypeScript merged with the map one and rejected as conflicting.
- The `@ratmath/complex` range is now `^0.3.0` (was `^0.2.0`, which under 0.x
  rejects complex 0.3.0, so the package could not be installed from npm).
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
