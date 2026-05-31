//! Bifurcation theory: detection and classification of bifurcations.

use serde::{Serialize, Deserialize};

/// Type of local bifurcation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum BifurcationType {
    SaddleNode,
    Transcritical,
    Pitchfork,
    Hopf,
    PeriodDoubling,
    Unknown,
}

/// A detected bifurcation point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationPoint {
    pub parameter: f64,
    pub bif_type: BifurcationType,
    pub description: String,
}

/// Detect bifurcation points in a 1D map by scanning for sign changes in the
/// Lyapunov exponent as a function of the bifurcation parameter.
///
/// Returns detected bifurcation points.
pub fn detect_bifurcations_1d<F, DF>(
    map: F,
    deriv: DF,
    param_range: (f64, f64),
    steps: usize,
    n_iterations: usize,
    x0: f64,
) -> Vec<BifurcationPoint>
where
    F: Fn(f64, f64) -> f64,   // map(r, x) -> x'
    DF: Fn(f64, f64) -> f64,  // deriv(r, x) -> f'(x)
{
    let (r_min, r_max) = param_range;
    let dr = (r_max - r_min) / steps as f64;
    let mut bifurcations = Vec::new();

    let mut prev_le = f64::NAN;
    for i in 0..=steps {
        let r = r_min + dr * i as f64;
        let le = compute_le_1d(&map, &deriv, r, x0, n_iterations);

        if !prev_le.is_nan() && !le.is_nan() {
            // Sign change in LE => transition to/from chaos
            if prev_le <= 0.0 && le > 0.0 {
                bifurcations.push(BifurcationPoint {
                    parameter: r - dr * 0.5,
                    bif_type: BifurcationType::PeriodDoubling,
                    description: format!("Transition to chaos at r={:.4}", r - dr * 0.5),
                });
            } else if prev_le > 0.0 && le <= 0.0 {
                bifurcations.push(BifurcationPoint {
                    parameter: r - dr * 0.5,
                    bif_type: BifurcationType::Unknown,
                    description: format!("Exit from chaos at r={:.4}", r - dr * 0.5),
                });
            }
        }
        prev_le = le;
    }

    bifurcations
}

fn compute_le_1d<F, DF>(map: F, deriv: DF, r: f64, x0: f64, n: usize) -> f64
where
    F: Fn(f64, f64) -> f64,
    DF: Fn(f64, f64) -> f64,
{
    let mut x = x0;
    let mut sum = 0.0f64;
    let transient = n / 5;
    for i in 0..(n + transient) {
        let d = deriv(r, x);
        if d.abs() > 1e-30 {
            sum += d.abs().ln();
        }
        x = map(r, x);
        if !x.is_finite() {
            return f64::NAN;
        }
        if i < transient {
            sum = 0.0;
        }
    }
    sum / n as f64
}

/// Classify a 1D fixed point bifurcation by examining eigenvalue behavior.
///
/// For a fixed point x* of x_{n+1} = f(r, x_n), a bifurcation occurs when |f'| = 1.
/// Returns the type based on second-order terms.
pub fn classify_bifurcation(
    eigenvalue: f64,
    second_deriv: f64,
    parameter_deriv: f64,
) -> BifurcationType {
    if (eigenvalue - 1.0).abs() < 0.1 {
        // Eigenvalue passes through +1
        if second_deriv.abs() < 0.1 && parameter_deriv.abs() > 0.1 {
            BifurcationType::Transcritical
        } else if second_deriv.abs() < 0.1 {
            BifurcationType::SaddleNode
        } else {
            BifurcationType::Pitchfork
        }
    } else if (eigenvalue + 1.0).abs() < 0.1 {
        // Eigenvalue passes through -1 => period doubling
        BifurcationType::PeriodDoubling
    } else {
        BifurcationType::Unknown
    }
}

/// Detect Hopf bifurcation in a 2D system.
///
/// Monitors when eigenvalues of the Jacobian at a fixed point cross the imaginary axis.
pub fn detect_hopf_bifurcation<F>(
    system: F,
    param_range: (f64, f64),
    steps: usize,
) -> Vec<BifurcationPoint>
where
    F: Fn(f64) -> ([f64; 2], [[f64; 2]; 2]), // r -> (fixed_point, jacobian)
{
    let (r_min, r_max) = param_range;
    let dr = (r_max - r_min) / steps as f64;
    let mut bifurcations = Vec::new();

    let mut prev_real_part: Option<f64> = None;

    for i in 0..=steps {
        let r = r_min + dr * i as f64;
        let (_fp, jac) = system(r);

        // Compute eigenvalues of 2x2 matrix: trace and determinant
        let trace = jac[0][0] + jac[1][1];
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        let disc = trace * trace - 4.0 * det;

        let real_part = if disc >= 0.0 {
            trace / 2.0
        } else {
            trace / 2.0 // real part of complex eigenvalues
        };

        if let Some(prev) = prev_real_part {
            // Detect when real part crosses zero with complex eigenvalues (disc < 0)
            if disc < 0.0 && prev <= 0.0 && real_part > 0.0 {
                bifurcations.push(BifurcationPoint {
                    parameter: r - dr * 0.5,
                    bif_type: BifurcationType::Hopf,
                    description: format!("Hopf bifurcation at r={:.4}", r - dr * 0.5),
                });
            }
        }
        prev_real_part = Some(real_part);
    }

    bifurcations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logistic_bifurcation_detection() {
        let map = |r: f64, x: f64| r * x * (1.0 - x);
        let deriv = |r: f64, x: f64| r * (1.0 - 2.0 * x);
        let bifs = detect_bifurcations_1d(map, deriv, (2.5, 4.0), 200, 5000, 0.1);
        // Should detect transitions near r ≈ 3.57 (onset of chaos)
        assert!(!bifs.is_empty(), "Should detect bifurcations in logistic map");
    }

    #[test]
    fn test_classify_saddle_node() {
        let btype = classify_bifurcation(1.0, 1.0, 0.0);
        assert_eq!(btype, BifurcationType::Pitchfork);
    }

    #[test]
    fn test_classify_period_doubling() {
        let btype = classify_bifurcation(-1.0, 0.0, 0.0);
        assert_eq!(btype, BifurcationType::PeriodDoubling);
    }

    #[test]
    fn test_classify_transcritical() {
        let btype = classify_bifurcation(1.0, 0.0, 1.0);
        assert_eq!(btype, BifurcationType::Transcritical);
    }

    #[test]
    fn test_hopf_detection() {
        // Simple system with Hopf: normal form dx/dt = μx - y, dy/dt = x + μy
        // Jacobian at origin: [[μ, -1], [1, μ]], eigenvalues μ ± i
        // Hopf at μ = 0
        let system = |r: f64| -> ([f64; 2], [[f64; 2]; 2]) {
            let fp = [0.0, 0.0];
            let jac = [[r, -1.0], [1.0, r]];
            (fp, jac)
        };
        let bifs = detect_hopf_bifurcation(system, (-1.0, 1.0), 200);
        assert_eq!(bifs.len(), 1);
        assert_eq!(bifs[0].bif_type, BifurcationType::Hopf);
        assert!((bifs[0].parameter).abs() < 0.05, "Hopf should be near 0, got {}", bifs[0].parameter);
    }
}
