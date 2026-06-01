# lau-chaos-theory

Chaos theory and nonlinear dynamics for Rust: Lyapunov exponents, fractal dimensions, strange attractors, bifurcation analysis, Poincaré sections, symbolic dynamics, period-doubling cascades, fractal generation, and agent unpredictability scoring.

## What This Does

This crate gives you the mathematical tools to analyze deterministic chaos. Feed it a dynamical system (a map or a flow), and it will tell you whether it's chaotic, measure how chaotic, characterize the attractor's geometry, detect bifurcation points, and produce fractal imagery. It also includes a higher-level "agent unpredictability" module that takes a raw time series and returns a scored report—useful for quantifying how unpredictable an agent's behavior is.

The library covers nine modules:

| Module | Purpose |
|---|---|
| `sensitivity` | Lyapunov exponents (1-D maps, n-D maps, continuous flows), Kaplan-Yorke dimension |
| `fractal_dimension` | Box-counting, correlation (Grassberger-Procaccia), information dimension |
| `bifurcation` | Bifurcation detection and classification (saddle-node, transcritical, pitchfork, Hopf, period-doubling) |
| `attractors` | Integrators for Lorenz, Rössler, Hénon-Heiles, Chua's circuit, and the Hénon map |
| `poincare` | Poincaré sections (directional and bidirectional) and return maps |
| `symbolic` | Partition encoding, shift spaces, topological entropy, Shannon entropy rate |
| `period_doubling` | Feigenbaum constants, logistic-map cascade analysis, superstable orbit finder |
| `fractal_gen` | Mandelbrot set, Julia sets, Burning Ship fractal (point tests and grids) |
| `agent_unpredictability` | Time-series Lyapunov estimation, behavioral regime classification, unpredictability scores |

**59 tests** across all modules.

## Key Idea

Chaos is deterministic but unpredictable. This crate operationalizes that paradox: Lyapunov exponents quantify exponential divergence of nearby trajectories; fractal dimensions measure the geometry of strange attractors; symbolic dynamics reduces continuous orbits to discrete sequences whose entropy bounds predictability. Every function is pure (no side effects, no global state), making it straightforward to compose analyses.

## Install

```toml
[dependencies]
lau-chaos-theory = "0.1.0"
```

Requires **Rust 2021 edition**. Dependencies:

- `nalgebra` 0.33 — linear algebra (matrices, QR decomposition)
- `serde` 1 (with `derive`) — serialization of results and parameters

## Quick Start

### Lyapunov exponent of the logistic map

```rust
use lau_chaos_theory::lyapunov_exponent_map;

let f  = |x: f64| 4.0 * x * (1.0 - x);   // logistic map, r = 4
let df = |x: f64| 4.0 * (1.0 - 2.0 * x); // derivative

let le = lyapunov_exponent_map(f, df, 0.1, 50_000);
println!("Lyapunov exponent = {:.4}  (theoretical: ln 2 ≈ 0.693)", le);
```

### Integrate the Lorenz attractor

```rust
use lau_chaos_theory::{lorenz_integrate, LorenzParams};

let params = LorenzParams::default(); // σ=10, ρ=28, β=8/3
let trajectory = lorenz_integrate(&params, [1.0, 1.0, 1.0], 0.01, 10_000);
println!("Final state: {:?}", trajectory.last().unwrap());
```

### Detect bifurcations in the logistic map

```rust
use lau_chaos_theory::detect_bifurcations_1d;

let map   = |r: f64, x: f64| r * x * (1.0 - x);
let deriv = |r: f64, x: f64| r * (1.0 - 2.0 * x);

let bifs = detect_bifurcations_1d(map, deriv, (2.5, 4.0), 200, 5000, 0.1);
for b in &bifs {
    println!("{:?} at r = {:.4}: {}", b.bif_type, b.parameter, b.description);
}
```

### Generate a Mandelbrot grid

```rust
use lau_chaos_theory::mandelbrot_grid;

let grid = mandelbrot_grid((-2.0, 1.0), (-1.5, 1.5), 800, 800, 256);
// grid[row][col] = Some(escape_iter) if escaped, None if in set
```

### Analyze agent unpredictability

```rust
use lau_chaos_theory::{analyze_scalar_unpredictability, build_agent_profile};

let series: Vec<f64> = /* your time series */;
let report = analyze_scalar_unpredictability(&series, vec![series_mean]);
println!("Chaotic? {}  Score: {:.3}", report.is_chaotic, report.unpredictability_score);

let profile = build_agent_profile(&series);
println!("Regime: {:?}", profile.behavioral_regime);
```

## API Reference

### Sensitivity (`sensitivity`)

| Function | Signature | Returns |
|---|---|---|
| `lyapunov_exponent_map` | `(f: Fn(f64)->f64, df: Fn(f64)->f64, x0, n) -> f64` | Lyapunov exponent for a 1-D map |
| `lyapunov_exponents_map` | `(f, _jacobian, x0: &DVector, n) -> Vec<f64>` | Full spectrum for an n-D map (QR method) |
| `lyapunov_exponents_flow` | `(f, x0: &DVector, dt, n_steps, renorm_interval) -> Vec<f64>` | Full spectrum for a continuous flow (RK4 + QR) |
| `max_lyapunov` | `(&[f64]) -> f64` | Largest exponent |
| `is_chaotic` | `(&[f64]) -> bool` | `true` if max exponent > 0 |
| `kaplan_yorke_dimension` | `(&[f64]) -> f64` | KY dimension estimate |

### Fractal Dimension (`fractal_dimension`)

| Function | Signature | Returns |
|---|---|---|
| `box_counting_dimension` | `(&[(f64,f64)], &[f64]) -> f64` | D_box via log-log regression |
| `correlation_dimension` | `(&[(f64,f64)], &[f64]) -> f64` | D_corr (Grassberger-Procaccia) |
| `information_dimension` | `(&[(f64,f64)], &[f64]) -> f64` | D_1 via Shannon entropy scaling |

### Bifurcation (`bifurcation`)

| Function | Signature | Returns |
|---|---|---|
| `detect_bifurcations_1d` | `(map, deriv, param_range, steps, n_iter, x0) -> Vec<BifurcationPoint>` | Scans for LE sign changes |
| `classify_bifurcation` | `(eigenvalue, second_deriv, param_deriv) -> BifurcationType` | Classifies by eigenvalue/curvature |
| `detect_hopf_bifurcation` | `(system, param_range, steps) -> Vec<BifurcationPoint>` | Detects Hopf in 2-D systems |

`BifurcationType`: `SaddleNode`, `Transcritical`, `Pitchfork`, `Hopf`, `PeriodDoubling`, `Unknown`.

### Attractors (`attractors`)

| Function | Params struct | Dimensions |
|---|---|---|
| `lorenz_integrate` | `LorenzParams { sigma, rho, beta }` | 3-D |
| `rossler_integrate` | `RosslerParams { a, b, c }` | 3-D |
| `henon_heiles_integrate` | `HenonHeilesParams { lambda }` | 4-D (x, y, px, py) |
| `chua_integrate` | `ChuaParams { alpha, beta, m0, m1 }` | 3-D |
| `henon_map` | (a, b inline) | 2-D discrete |

All integrators use RK4. Params implement `Default` with canonical chaotic values.

### Poincaré Sections (`poincare`)

| Function | Description |
|---|---|
| `poincare_section(trajectory, section_fn, offset)` | All plane crossings with interpolation |
| `poincare_section_directed(trajectory, section_fn, offset, positive_direction)` | Direction-filtered crossings |
| `return_map(crossings, component)` | Successive-intersection scatter pairs |

`const N: usize` generics allow any dimension.

### Symbolic Dynamics (`symbolic`)

| Item | Description |
|---|---|
| `Partition::new(boundaries)` | Creates a state-space partition |
| `encode_itinerary(trajectory, partition)` | Maps scalar series → symbol sequence |
| `ShiftSpace { n_symbols, forbidden }` | Defines allowed/forbidden transitions |
| `topological_entropy(sequence, block_size)` | Estimates h_top via n-block counting |
| `shannon_entropy_rate(sequence)` | Transition-entropy rate |

### Period Doubling (`period_doubling`)

| Item | Description |
|---|---|
| `FEIGENBAUM_DELTA` (4.6692…) | First Feigenbaum constant |
| `FEIGENBAUM_ALPHA` (2.5029…) | Second Feigenbaum constant |
| `logistic_bifurcation_points(r_min, r_max, resolution)` | Finds cascade bifurcation values |
| `estimate_feigenbaum_delta(&[f64])` | Computes δ from bifurcation-point list |
| `verify_feigenbaum_delta()` | Verifies δ from known logistic points |
| `superstable_parameter(period_power)` | Finds r where 2ⁿ-orbit is superstable |

### Fractal Generation (`fractal_gen`)

| Function | Description |
|---|---|
| `mandelbrot(c_re, c_im, max_iter)` | Point test → `FractalResult` |
| `julia(z_re, z_im, c_re, c_im, max_iter)` | Point test → `FractalResult` |
| `burning_ship(c_re, c_im, max_iter)` | Point test → `FractalResult` |
| `mandelbrot_grid(x_range, y_range, w, h, max_iter)` | Grid of escape times |
| `julia_grid(c_re, c_im, x_range, y_range, w, h, max_iter)` | Grid of escape times |

`FractalResult { in_set, iterations, smooth_escape }` — smooth escape uses the normalized iteration count for band-free coloring.

### Agent Unpredictability (`agent_unpredictability`)

| Function | Returns |
|---|---|
| `analyze_scalar_unpredictability(series, boundaries)` | `UnpredictabilityReport` with LE, entropy rates, score |
| `build_agent_profile(series)` | `AgentProfile` with variance, ACF(1), regime classification |

`BehavioralRegime`: `FixedPoint`, `Periodic`, `Quasiperiodic`, `Chaotic`, `Unknown`.

## How It Works

### Lyapunov Exponents

For 1-D maps, the Lyapunov exponent is computed directly from the definition:

$$\lambda = \frac{1}{N} \sum_{n=0}^{N-1} \ln|f'(x_n)|$$

For n-dimensional systems, the algorithm maintains a QR decomposition of the evolving tangent-space basis. At each step, the Jacobian multiplies the current orthonormal frame, and periodic QR re-orthonormalization prevents collapse. The logarithms of the diagonal elements of R accumulate the Lyapunov sums.

Continuous flows use RK4 integration of both the state and the variational equations simultaneously, with QR renormalization at configurable intervals.

### Fractal Dimensions

**Box-counting** divides the plane into ε-boxes, counts occupied boxes N(ε) at multiple scales, and fits a line to ln N(ε) vs ln(1/ε). The slope is D_box.

**Correlation dimension** (Grassberger-Procaccia) counts pairs within distance ε, forming the correlation integral C(ε), and takes the slope of ln C(ε) vs ln ε.

**Information dimension** computes the Shannon entropy H(ε) of the box occupation distribution at each scale and takes the slope of H(ε) vs ln(1/ε).

### Bifurcation Detection

The library scans a parameter range and computes the Lyapunov exponent at each value. Sign changes (negative → positive or vice versa) mark transitions to/from chaos. Classification uses eigenvalue position: crossing +1 gives saddle-node/transcritical/pitchfork (distinguished by second derivative and parameter derivative), crossing −1 gives period-doubling. Hopf detection monitors the real part of complex eigenvalue pairs in 2-D Jacobians.

### Attractor Integration

All continuous systems use a hand-rolled RK4 integrator. This is intentional—the crate avoids depending on ODE solver crates to stay self-contained. The timestep and step count are caller-controlled, letting you trade speed for accuracy.

### Poincaré Sections

Given a trajectory and a hyperplane (defined by a normal vector and offset), the algorithm detects sign changes in the dot product `⟨state, normal⟩ − offset` between consecutive points. Linear interpolation recovers the crossing point to sub-step accuracy. Direction-aware variants filter by sign of the crossing.

### Symbolic Dynamics

A `Partition` divides the real line into regions. `encode_itinerary` maps a scalar trajectory to a symbol sequence. Topological entropy is estimated by counting distinct n-blocks (finite-data approximation). Shannon entropy rate comes from the transition probability matrix of the symbol sequence.

### Period Doubling and Feigenbaum Constants

The logistic map's period-doubling cascade is detected by iterating at each parameter value and checking orbit periodicity. The Feigenbaum δ constant (ratio of successive bifurcation intervals) can be estimated from detected points and compared to the universal value 4.6692…

### Fractal Generation

Uses the escape-time algorithm: iterate z → z² + c (or variants) until |z| > 2 or max iterations reached. Smooth coloring uses the normalized iteration count:

$$\mu = n + 1 - \frac{\ln \ln |z_n|}{\ln 2}$$

The Burning Ship variant takes absolute values of real and imaginary parts before squaring.

### Agent Unpredictability

This module applies the lower-level tools to raw time series. It estimates a Lyapunov exponent using a nearest-neighbor divergence method (Rosenstein et al.), computes Shannon and topological entropy on a symbolic encoding, and combines them into a single unpredictability score. `build_agent_profile` adds variance, autocorrelation, and regime classification.

## The Math

### Lyapunov Exponents

The maximal Lyapunov exponent (MLE) measures the average exponential rate of divergence of nearby trajectories:

$$\lambda_{\max} = \lim_{t \to \infty} \frac{1}{t} \ln \frac{|\delta(t)|}{|\delta(0)|}$$

- **λ > 0**: chaos (exponential sensitivity)
- **λ = 0**: marginal (periodic or quasiperiodic)
- **λ < 0**: stable (convergence to attractor)

For the logistic map at r = 4: λ = ln 2 ≈ 0.693. For the Lorenz system with standard parameters: λ₁ ≈ 0.906, λ₂ = 0, λ₃ ≈ −14.57.

### Kaplan-Yorke Dimension

Given the Lyapunov spectrum sorted in decreasing order λ₁ ≥ λ₂ ≥ … ≥ λ_n, the KY dimension is:

$$D_{KY} = j + \frac{\sum_{i=1}^{j} \lambda_i}{|\lambda_{j+1}|}$$

where j is the largest index with cumulative sum still positive. For Lorenz: D_KY ≈ 2 + 0.906/14.57 ≈ 2.062.

### Feigenbaum Universality

The period-doubling cascade of unimodal maps converges geometrically with ratio:

$$\delta = \lim_{n \to \infty} \frac{r_n - r_{n-1}}{r_{n+1} - r_n} = 4.6692\ldots$$

This constant is universal—it's the same for any map with a quadratic maximum, not just the logistic map.

### Fractal Dimensions

For a self-similar set, the box-counting dimension satisfies:

$$N(\varepsilon) \sim \varepsilon^{-D}$$

The correlation dimension generalizes this by considering pair distances:

$$C(\varepsilon) = \frac{2}{N(N-1)} \sum_{i < j} \Theta(\varepsilon - \|x_i - x_j\|)$$

In general: D₀ (box) ≤ D₁ (information) ≤ D₂ (correlation), with equality for homogeneous sets.

### Topological Entropy

For a shift space with forbidden transitions described by a transition matrix A, the topological entropy equals:

$$h_{\text{top}} = \ln \rho(A)$$

where ρ(A) is the spectral radius. For a full n-shift: h_top = ln n.

## License

MIT
