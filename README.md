# @ratmath/chaos

Dynamical systems and chaos for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/). Docs: [ratmath.com/chaos](https://ratmath.com/chaos/).

A chaotic map doubles every rounding error at each step, so after about fifty
steps a floating-point orbit has no correct digit left and nothing warns you.
Here an orbit is the true one, exactly or inside an interval; periodic points
are algebraic numbers; and the parameters at which cycles are born or double
are proven.

```js
import { logistic, lyapunov, bifurcationPoint, itinerary } from "@ratmath/chaos";

const f = logistic(4);                       // x → 4x(1 − x)
f.orbit("1/3", 3);                           // 1/3, 8/9, 32/81, 6272/6561: the true orbit
f.compare("1/10", 80).parted;                // 52: the step where doubles have no digit right
f.periodicPoints(2).map((p) => p.text);      // ["(5 − √5)/8", "(5 + √5)/8"]
lyapunov(f, "1/10", { steps: 300 }).value;   // an interval, along the true orbit
bifurcationPoint(logistic(), { parameter: "r", period: 2, near: [3.44, 0.44] }).parameter;  // 1 + √6, proven
itinerary(f, "1/10", 12);                    // "LLRLRRRLLRLR"
```

## What is built

| Area | Functions |
|---|---|
| Maps and orbits | `map(formula, { params })`, `logistic`, `tent`, `quadratic`, `doubling`; `.orbit` (exact), `.enclosedOrbit` (certified), `.floatOrbit`, `.compare`, `.cobweb` |
| Periodic points | `.fixedPoints()`, `.periodicPoints(n)`, `.cycles(n)` as exact algebraic numbers with multiplier and stability; `.periodicPointsIn(n, [a, b])` certified for any formula |
| Exponents | `lyapunov` (certified average along the true orbit), `sensitivity` |
| Bifurcations | `bifurcationPoint`, `periodDoubling`, `periodicWindow`, `superstable` (all proven); `bifurcationDiagram` (estimated) |
| Symbolic dynamics | `itinerary`, `kneading`, `sharkovskii` |
| Maps of the plane | `planeMap`, `henon`, `standardMap`: exact, enclosed and float orbits; `.periodicPoints(n, region)` each proven and the list proven complete; `.exactPeriodicPoints(n)` as algebraic numbers; `catMap`, `torusMap` (exact periods, counts, entropy), `horseshoe` (symbolic dynamics in exact coordinates) |
| Flows | `flow(field)`, `lorenz`, `rossler`, `vanDerPol`, `pendulum`: `.solve` encloses the true solution; `.returnMap` and `.periodicOrbit` prove closed orbits with period and multipliers; `.equilibria(region)` finds and classifies every rest point; `.bifurcations` brackets saddle-node, pitchfork, transcritical and Hopf parameters |
| Fractals | `IFS`, `Affine`, `attractors` (exact stages, contraction bounds, similarity dimension exact or certified, chaos game estimated), `moran`, `boxCounting`; `cantor({ base, digits })` with exact membership of fractions |
| L-systems | `lsystem`, `lsystems`: symbol counts by matrix power; turtle paths exact in a + b√3 at multiples of 30° |
| Cellular automata | `elementary(rule)` (equivalents, additivity, cycles), `life(cells, { rule })` with `classify` (still life, oscillator, spaceship and speed) |
| Escape-time pictures | `mandelbrot`, `julia`, `newtonBasins`: every pixel proved inside, proved outside or undecided; `mandelbrotPoint`, `juliaPoint` for exact points |

Three kinds of answer, each labelled: exact (rationals and algebraic numbers),
certified (intervals proven to contain the value) and estimated (pictures and
observed limits).

Every phase of the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#10-a-dynamical-systems-package-ratmathchaos)
for this package is built; the pictures are in `@ratmath/graphics`.

## Layout

- `index.js`, `src/`: the library (`tree.js` evaluators, `map.js` maps and
  periodic points, `dynamics.js` exponents, bifurcations and symbols,
  `fractals.js`, `lsystem.js`, `automata.js`, `escape.js`, `plane.js` maps of
  the plane, `fi.js` outward-rounded float intervals, Taylor series and
  Krawczyk tests, `flow.js` the enclosing solver, sections, periodic orbits,
  equilibria and bifurcations).
- `vendor/analytic`, `vendor/algebraic`, `vendor/complex`: git submodules, used as workspaces
  until those packages are on npm. Clone with `--recurse-submodules`.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages; `bun run docs:build` writes `docs/*.html`,
  typesets the formulas and runs every page's code box.

## License

MIT
