# lau-chaos-theory

> Chaos theory and nonlinear dynamics: Lyapunov exponents, fractals, strange attractors, bifurcation analysis, symbolic dynamics, and agent unpredictability

## What This Does

Chaos theory and nonlinear dynamics: Lyapunov exponents, fractals, strange attractors, bifurcation analysis, symbolic dynamics, and agent unpredictability. Part of the PLATO/LAU ecosystem — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## The Key Idea

This crate implements the core abstractions needed for its domain, with a focus on correctness, composability, and conservation guarantees. Every public type is serializable (serde), every algorithm is tested, and every invariant is verified.

## Install

```bash
cargo add lau-chaos-theory
```

## Quick Start

See the API Reference below for complete usage. Key entry points:

```rust
use lau_chaos_theory::*;
// See types and methods below for complete usage
```

## API Reference

```rust
pub struct LorenzParams 
pub fn lorenz_integrate(
pub struct RosslerParams 
pub fn rossler_integrate(
pub struct HenonHeilesParams 
pub fn henon_heiles_integrate(
pub struct ChuaParams 
pub fn chua_integrate(
pub fn henon_map(a: f64, b: f64, x0: [f64; 2], n: usize) -> Vec<[f64; 2]> 
pub enum BifurcationType 
pub struct BifurcationPoint 
pub fn detect_bifurcations_1d<F, DF>(
pub fn classify_bifurcation(
pub fn detect_hopf_bifurcation<F>(
pub fn poincare_section<const N: usize>(
pub fn poincare_section_directed<const N: usize>(
pub fn return_map(crossings: &[[f64; 2]], component: usize) -> Vec<(f64, f64)> 
pub struct FractalResult 
pub fn mandelbrot(c_re: f64, c_im: f64, max_iter: usize) -> FractalResult 
pub fn julia(z_re: f64, z_im: f64, c_re: f64, c_im: f64, max_iter: usize) -> FractalResult 
pub fn mandelbrot_grid(
pub fn julia_grid(
pub fn burning_ship(c_re: f64, c_im: f64, max_iter: usize) -> FractalResult 
pub struct Partition 
    pub fn new(boundaries: Vec<f64>) -> Self 
    pub fn symbol(&self, x: f64) -> usize 
    pub fn n_symbols(&self) -> usize 
pub fn encode_itinerary(trajectory: &[f64], partition: &Partition) -> Vec<usize> 
pub struct ShiftSpace 
    pub fn full(n: usize) -> Self 
    pub fn is_allowed(&self, from: usize, to: usize) -> bool 
    pub fn allowed_transitions(&self, from: usize) -> usize 
pub fn topological_entropy(sequence: &[usize], block_size: usize) -> f64 
pub fn shannon_entropy_rate(sequence: &[usize]) -> f64 
pub enum BehaviorMetric 
pub struct UnpredictabilityReport 
pub fn analyze_scalar_unpredictability(
pub struct AgentProfile 
pub enum BehavioralRegime 
pub fn build_agent_profile(series: &[f64]) -> AgentProfile 
pub fn lyapunov_exponent_map<F, G>(f: F, df: G, x0: f64, n: usize) -> f64
pub fn lyapunov_exponents_map<F, G>(
pub fn lyapunov_exponents_flow<F>(
pub fn max_lyapunov(exponents: &[f64]) -> f64 
pub fn is_chaotic(exponents: &[f64]) -> bool 
pub fn kaplan_yorke_dimension(exponents: &[f64]) -> f64 
pub fn box_counting_dimension(points: &[(f64, f64)], scales: &[f64]) -> f64 
pub fn correlation_dimension(points: &[(f64, f64)], scales: &[f64]) -> f64 
pub fn information_dimension(points: &[(f64, f64)], scales: &[f64]) -> f64 
pub fn logistic_bifurcation_points(r_min: f64, r_max: f64, resolution: usize) -> Vec<f64> 
pub fn estimate_feigenbaum_delta(bif_points: &[f64]) -> f64 
pub fn verify_feigenbaum_delta() -> f64 
pub fn superstable_parameter(period_power: u32) -> f64 
```

## How It Works

Read the source in `src/` for full implementation details. All algorithms are documented with inline comments explaining the mathematical foundations.

## The Math

This crate implements formal mathematical constructs. See the source documentation for theorem statements and proofs of correctness.

## Testing

**59 tests** covering construction, serialization, correctness properties, edge cases, and composability with other lau-* crates.

## License

MIT
