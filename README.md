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

Three kinds of answer, each labelled: exact (rationals and algebraic numbers),
certified (intervals proven to contain the value) and estimated (pictures and
observed limits).

## Still planned

| Phase | Contents |
|---|---|
| 3 | Maps of the plane; flows, sections, periodic orbits, bifurcations |
| 4 | Iterated function systems, Cantor sets, dimensions, L-systems, cellular automata |
| 5 | Julia and Mandelbrot sets, Newton basins |
| 6 | Pictures |

The plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#10-a-dynamical-systems-package-ratmathchaos),
and `plan` exports the list.

## Layout

- `index.js`, `src/`: the library (`tree.js` evaluators, `map.js` maps and
  periodic points, `dynamics.js` exponents, bifurcations and symbols).
- `vendor/analytic`, `vendor/algebraic`: git submodules, used as workspaces
  until those packages are on npm. Clone with `--recurse-submodules`.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages; `bun run docs:build` writes `docs/*.html`,
  typesets the formulas and runs every page's code box.

## License

MIT
