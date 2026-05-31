//! Strange attractors: Lorenz, Rössler, Hénon-Heiles, Chua's circuit.

use serde::{Serialize, Deserialize};

/// Lorenz system parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LorenzParams {
    pub sigma: f64,
    pub rho: f64,
    pub beta: f64,
}

impl Default for LorenzParams {
    fn default() -> Self {
        Self { sigma: 10.0, rho: 28.0, beta: 8.0 / 3.0 }
    }
}

/// Integrate the Lorenz system using RK4.
pub fn lorenz_integrate(
    params: &LorenzParams,
    x0: [f64; 3],
    dt: f64,
    steps: usize,
) -> Vec<[f64; 3]> {
    let LorenzParams { sigma, rho, beta } = *params;
    let mut state = x0;
    let mut trajectory = Vec::with_capacity(steps);

    for _ in 0..steps {
        let k1 = lorenz_rhs(sigma, rho, beta, &state);
        let mut s2 = state;
        for i in 0..3 { s2[i] += dt * 0.5 * k1[i]; }
        let k2 = lorenz_rhs(sigma, rho, beta, &s2);
        let mut s3 = state;
        for i in 0..3 { s3[i] += dt * 0.5 * k2[i]; }
        let k3 = lorenz_rhs(sigma, rho, beta, &s3);
        let mut s4 = state;
        for i in 0..3 { s4[i] += dt * k3[i]; }
        let k4 = lorenz_rhs(sigma, rho, beta, &s4);

        for i in 0..3 {
            state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        trajectory.push(state);
    }

    trajectory
}

fn lorenz_rhs(sigma: f64, rho: f64, beta: f64, s: &[f64; 3]) -> [f64; 3] {
    [
        sigma * (s[1] - s[0]),
        s[0] * (rho - s[2]) - s[1],
        s[0] * s[1] - beta * s[2],
    ]
}

/// Rössler system parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RosslerParams {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl Default for RosslerParams {
    fn default() -> Self {
        Self { a: 0.2, b: 0.2, c: 5.7 }
    }
}

/// Integrate the Rössler system using RK4.
pub fn rossler_integrate(
    params: &RosslerParams,
    x0: [f64; 3],
    dt: f64,
    steps: usize,
) -> Vec<[f64; 3]> {
    let RosslerParams { a, b, c } = *params;
    let mut state = x0;
    let mut trajectory = Vec::with_capacity(steps);

    for _ in 0..steps {
        let k1 = rossler_rhs(a, b, c, &state);
        let mut s2 = state;
        for i in 0..3 { s2[i] += dt * 0.5 * k1[i]; }
        let k2 = rossler_rhs(a, b, c, &s2);
        let mut s3 = state;
        for i in 0..3 { s3[i] += dt * 0.5 * k2[i]; }
        let k3 = rossler_rhs(a, b, c, &s3);
        let mut s4 = state;
        for i in 0..3 { s4[i] += dt * k3[i]; }
        let k4 = rossler_rhs(a, b, c, &s4);

        for i in 0..3 {
            state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        trajectory.push(state);
    }

    trajectory
}

fn rossler_rhs(a: f64, b: f64, c: f64, s: &[f64; 3]) -> [f64; 3] {
    [
        -s[1] - s[2],
        s[0] + a * s[1],
        b + s[2] * (s[0] - c),
    ]
}

/// Hénon-Heiles system (Hamiltonian system with chaotic behavior).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HenonHeilesParams {
    pub lambda: f64,
}

impl Default for HenonHeilesParams {
    fn default() -> Self {
        Self { lambda: 1.0 }
    }
}

/// Integrate Hénon-Heiles using RK4.
/// State: [x, y, px, py]
pub fn henon_heiles_integrate(
    params: &HenonHeilesParams,
    x0: [f64; 4],
    dt: f64,
    steps: usize,
) -> Vec<[f64; 4]> {
    let lambda = params.lambda;
    let mut state = x0;
    let mut trajectory = Vec::with_capacity(steps);

    for _ in 0..steps {
        let k1 = hh_rhs(lambda, &state);
        let mut s2 = state;
        for i in 0..4 { s2[i] += dt * 0.5 * k1[i]; }
        let k2 = hh_rhs(lambda, &s2);
        let mut s3 = state;
        for i in 0..4 { s3[i] += dt * 0.5 * k2[i]; }
        let k3 = hh_rhs(lambda, &s3);
        let mut s4 = state;
        for i in 0..4 { s4[i] += dt * k3[i]; }
        let k4 = hh_rhs(lambda, &s4);

        for i in 0..4 {
            state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        trajectory.push(state);
    }

    trajectory
}

fn hh_rhs(lambda: f64, s: &[f64; 4]) -> [f64; 4] {
    // x, y, px, py
    let x = s[0]; let y = s[1]; let px = s[2]; let py = s[3];
    [
        px,
        py,
        -x - 2.0 * lambda * x * y,
        -y - lambda * (x * x - y * y),
    ]
}

/// Chua's circuit parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChuaParams {
    pub alpha: f64,
    pub beta: f64,
    pub m0: f64,
    pub m1: f64,
}

impl Default for ChuaParams {
    fn default() -> Self {
        Self { alpha: 15.6, beta: 28.0, m0: -1.143, m1: -0.714 }
    }
}

/// Chua's circuit nonlinearity.
fn chua_nonlinearity(x: f64, m0: f64, m1: f64) -> f64 {
    // Piecewise linear: m1*x + 0.5*(m0 - m1)*(|x+1| - |x-1|)
    m1 * x + 0.5 * (m0 - m1) * ((x + 1.0).abs() - (x - 1.0).abs())
}

/// Integrate Chua's circuit using RK4.
pub fn chua_integrate(
    params: &ChuaParams,
    x0: [f64; 3],
    dt: f64,
    steps: usize,
) -> Vec<[f64; 3]> {
    let ChuaParams { alpha, beta, m0, m1 } = *params;
    let mut state = x0;
    let mut trajectory = Vec::with_capacity(steps);

    for _ in 0..steps {
        let k1 = chua_rhs(alpha, beta, m0, m1, &state);
        let mut s2 = state;
        for i in 0..3 { s2[i] += dt * 0.5 * k1[i]; }
        let k2 = chua_rhs(alpha, beta, m0, m1, &s2);
        let mut s3 = state;
        for i in 0..3 { s3[i] += dt * 0.5 * k2[i]; }
        let k3 = chua_rhs(alpha, beta, m0, m1, &s3);
        let mut s4 = state;
        for i in 0..3 { s4[i] += dt * k3[i]; }
        let k4 = chua_rhs(alpha, beta, m0, m1, &s4);

        for i in 0..3 {
            state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        trajectory.push(state);
    }

    trajectory
}

fn chua_rhs(alpha: f64, beta: f64, m0: f64, m1: f64, s: &[f64; 3]) -> [f64; 3] {
    let h = chua_nonlinearity(s[0], m0, m1);
    [
        alpha * (s[1] - s[0] - h),
        s[0] - s[1] + s[2],
        -beta * s[1],
    ]
}

/// Hénon map (discrete 2D attractor).
pub fn henon_map(a: f64, b: f64, x0: [f64; 2], n: usize) -> Vec<[f64; 2]> {
    let mut state = x0;
    let mut trajectory = Vec::with_capacity(n);
    for _ in 0..n {
        let xn = 1.0 - a * state[0] * state[0] + b * state[1];
        let yn = state[0];
        state = [xn, yn];
        if !state[0].is_finite() { break; }
        trajectory.push(state);
    }
    trajectory
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lorenz_bounded() {
        let params = LorenzParams::default();
        let traj = lorenz_integrate(&params, [1.0, 1.0, 1.0], 0.01, 10000);
        // Lorenz attractor should be bounded
        for s in &traj {
            assert!(s[0].abs() < 50.0, "x should be bounded, got {}", s[0]);
            assert!(s[1].abs() < 50.0, "y should be bounded, got {}", s[1]);
            assert!(s[2] < 60.0 && s[2] > -5.0, "z should be bounded, got {}", s[2]);
        }
    }

    #[test]
    fn test_lorenz_sensitivity() {
        let params = LorenzParams::default();
        let t1 = lorenz_integrate(&params, [1.0, 1.0, 1.0], 0.01, 5000);
        let t2 = lorenz_integrate(&params, [1.0 + 1e-10, 1.0, 1.0], 0.01, 5000);
        // After some time, trajectories should diverge
        let d_start = ((t1[100][0] - t2[100][0]).powi(2) +
                       (t1[100][1] - t2[100][1]).powi(2) +
                       (t1[100][2] - t2[100][2]).powi(2)).sqrt();
        let d_end = ((t1[4999][0] - t2[4999][0]).powi(2) +
                     (t1[4999][1] - t2[4999][1]).powi(2) +
                     (t1[4999][2] - t2[4999][2]).powi(2)).sqrt();
        assert!(d_end > d_start * 100.0, "Should show sensitivity: d_start={:.6}, d_end={:.6}", d_start, d_end);
    }

    #[test]
    fn test_rossler_bounded() {
        let params = RosslerParams::default();
        let traj = rossler_integrate(&params, [0.1, 0.1, 0.1], 0.01, 50000);
        // Rössler should stay bounded (though with spikes)
        let max_z = traj.iter().map(|s| s[2]).fold(f64::NEG_INFINITY, f64::max);
        assert!(max_z < 50.0, "Rössler z should be bounded, max={}", max_z);
    }

    #[test]
    fn test_henon_map_attractor() {
        let traj = henon_map(1.4, 0.3, [0.0, 0.0], 5000);
        // Discard transient
        let pts: Vec<_> = traj.iter().skip(1000).collect();
        // Hénon attractor should be bounded
        for p in &pts {
            assert!(p[0].abs() < 5.0, "Hénon x should be bounded");
            assert!(p[1].abs() < 5.0, "Hénon y should be bounded");
        }
        assert!(pts.len() > 1000);
    }

    #[test]
    fn test_chua_double_scroll() {
        let params = ChuaParams::default();
        let traj = chua_integrate(&params, [0.1, 0.0, 0.1], 0.001, 50000);
        // Chua's circuit should exhibit double scroll
        // Both positive and negative x values should appear (after transient)
        let pts: Vec<_> = traj.iter().skip(10000).collect();
        let has_positive = pts.iter().any(|p| p[0] > 1.0);
        let has_negative = pts.iter().any(|p| p[0] < -1.0);
        assert!(has_positive && has_negative, "Chua should visit both scrolls");
    }

    #[test]
    fn test_henon_heiles_energy_conservation() {
        let params = HenonHeilesParams::default();
        let traj = henon_heiles_integrate(&params, [0.0, 0.0, 0.3, 0.3], 0.01, 10000);
        // Hénon-Heiles is Hamiltonian, energy should be approximately conserved
        let e0 = hh_energy(1.0, &traj[0]);
        let ef = hh_energy(1.0, &traj[9999]);
        assert!((e0 - ef).abs() / e0.abs() < 0.01, "Energy should be conserved: E0={}, Ef={}", e0, ef);
    }

    fn hh_energy(lambda: f64, s: &[f64; 4]) -> f64 {
        let x = s[0]; let y = s[1]; let px = s[2]; let py = s[3];
        0.5 * (px * px + py * py) + 0.5 * (x * x + y * y)
            + lambda * (x * x * y - y * y * y / 3.0)
    }
}
