//! Poincaré sections and return maps.

/// Compute a Poincaré section by detecting when a trajectory crosses a plane.
///
/// Detects zero crossings of `section_fn(state)` (sign change from + to - or vice versa)
/// and interpolates the crossing point.
pub fn poincare_section<const N: usize>(
    trajectory: &[[f64; N]],
    section_fn: &[f64; N],
    offset: f64,
) -> Vec<[f64; N]> {
    let mut crossings = Vec::new();

    for i in 1..trajectory.len() {
        let v0 = dot(&trajectory[i - 1], section_fn) - offset;
        let v1 = dot(&trajectory[i], section_fn) - offset;

        if v0 * v1 < 0.0 {
            // Linear interpolation to find crossing
            let t = v0 / (v0 - v1);
            let mut crossing = [0.0f64; N];
            for j in 0..N {
                crossing[j] = trajectory[i - 1][j] + t * (trajectory[i][j] - trajectory[i - 1][j]);
            }
            crossings.push(crossing);
        }
    }

    crossings
}

/// Direction-aware Poincaré section: only count crossings in a specific direction.
pub fn poincare_section_directed<const N: usize>(
    trajectory: &[[f64; N]],
    section_fn: &[f64; N],
    offset: f64,
    positive_direction: bool,
) -> Vec<[f64; N]> {
    let mut crossings = Vec::new();

    for i in 1..trajectory.len() {
        let v0 = dot(&trajectory[i - 1], section_fn) - offset;
        let v1 = dot(&trajectory[i], section_fn) - offset;

        let crosses = if positive_direction {
            v0 < 0.0 && v1 >= 0.0
        } else {
            v0 >= 0.0 && v1 < 0.0
        };

        if crosses {
            let t = if (v0 - v1).abs() > 1e-15 {
                v0 / (v0 - v1)
            } else {
                0.5
            };
            let mut crossing = [0.0f64; N];
            for j in 0..N {
                crossing[j] = trajectory[i - 1][j] + t * (trajectory[i][j] - trajectory[i - 1][j]);
            }
            crossings.push(crossing);
        }
    }

    crossings
}

/// Build a return map from Poincaré section data: maps successive intersection values.
pub fn return_map(crossings: &[[f64; 2]], component: usize) -> Vec<(f64, f64)> {
    let vals: Vec<f64> = crossings.iter().map(|c| c[component]).collect();
    let mut map = Vec::new();
    for i in 1..vals.len() {
        map.push((vals[i - 1], vals[i]));
    }
    map
}

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    let mut s = 0.0;
    for i in 0..N { s += a[i] * b[i]; }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::attractors::{lorenz_integrate, LorenzParams};

    #[test]
    fn test_poincare_lorenz_section() {
        let params = LorenzParams::default();
        let traj = lorenz_integrate(&params, [1.0, 1.0, 1.0], 0.005, 100000);
        // Section z = rho - 1 = 27
        let section_fn = [0.0, 0.0, 1.0];
        let crossings = poincare_section(&traj, &section_fn, 27.0);
        // Should have multiple crossings
        assert!(crossings.len() > 10, "Should have crossings, got {}", crossings.len());
        // All crossings should be near z=27
        for c in &crossings {
            assert!((c[2] - 27.0).abs() < 0.1, "Crossing z should be ~27, got {}", c[2]);
        }
    }

    #[test]
    fn test_poincare_directed() {
        let params = LorenzParams::default();
        let traj = lorenz_integrate(&params, [1.0, 1.0, 1.0], 0.005, 100000);
        let section_fn = [0.0, 0.0, 1.0];
        let crossings_up = poincare_section_directed(&traj, &section_fn, 27.0, true);
        let crossings_down = poincare_section_directed(&traj, &section_fn, 27.0, false);
        // Total up + down should equal total crossings
        let total = poincare_section(&traj, &section_fn, 27.0);
        assert_eq!(crossings_up.len() + crossings_down.len(), total.len());
    }

    #[test]
    fn test_return_map() {
        let crossings = vec![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]];
        let map = return_map(&crossings, 0);
        assert_eq!(map, vec![(1.0, 3.0), (3.0, 5.0)]);
    }

    #[test]
    fn test_poincare_simple_periodic() {
        // Simple circular orbit: x = cos(t), y = sin(t)
        let n = 10000;
        let dt = 0.001;
        let traj: Vec<[f64; 2]> = (0..n)
            .map(|i| {
                let t = i as f64 * dt;
                [t.cos(), t.sin()]
            })
            .collect();
        // Section y = 0, x > 0 (crossing upward through positive x-axis)
        let section_fn = [0.0, 1.0];
        let crossings = poincare_section_directed(&traj, &section_fn, 0.0, true);
        // Should get roughly n*dt/(2π) ≈ 1.59 crossings
        assert!(crossings.len() >= 1, "Should detect periodic crossing");
        for c in &crossings {
            assert!((c[0] - 1.0).abs() < 0.05, "x should be ~1 at y=0 crossing, got {}", c[0]);
        }
    }
}
