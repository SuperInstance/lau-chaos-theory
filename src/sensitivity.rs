//! Sensitivity analysis — Lyapunov exponents for maps and continuous systems.

use nalgebra::{DVector, DMatrix};

/// Compute the Lyapunov exponent for a 1D map x_{n+1} = f(x_n).
///
/// Uses the formula: λ = (1/N) Σ ln|f'(x_n)|
pub fn lyapunov_exponent_map<F, G>(f: F, df: G, x0: f64, n: usize) -> f64
where
    F: Fn(f64) -> f64,
    G: Fn(f64) -> f64,
{
    let mut x = x0;
    let mut sum = 0.0;
    let transient = n / 10;
    for i in 0..(n + transient) {
        let deriv = df(x);
        if deriv.abs() > 0.0 {
            sum += deriv.abs().ln();
        }
        x = f(x);
        if i >= transient && !x.is_finite() {
            return sum / (i - transient + 1) as f64;
        }
    }
    sum / n as f64
}

/// Compute all Lyapunov exponents for an n-dimensional map.
///
/// Uses QR decomposition method for the Jacobian product.
pub fn lyapunov_exponents_map<F, G>(
    f: F,
    _jacobian: G,
    x0: &DVector<f64>,
    n: usize,
) -> Vec<f64>
where
    F: Fn(&DVector<f64>) -> (DVector<f64>, DMatrix<f64>),
    G: Fn(&DVector<f64>) -> (DVector<f64>, DMatrix<f64>),
{
    let dim = x0.nrows();
    let mut x = x0.clone();
    let mut q = DMatrix::<f64>::identity(dim, dim);
    let mut sums = vec![0.0f64; dim];

    let transient = n / 10;

    for i in 0..(n + transient) {
        let (x_new, jac) = f(&x);
        x = x_new;
        let j_q = jac * &q;
        let qr = j_q.qr();
        let q_new = qr.q();
        let r = qr.r();
        for j in 0..dim {
            let diag = r[(j, j)].abs();
            if diag > 0.0 {
                sums[j] += diag.ln();
            }
        }
        q = q_new;

        if i >= transient && !x.iter().all(|v| v.is_finite()) {
            let count = (i - transient + 1) as f64;
            return sums.iter().map(|s| s / count).collect();
        }
    }

    sums.iter().map(|s| s / n as f64).collect()
}

/// Compute Lyapunov exponents for a continuous system (flow) using time integration.
///
/// `f` returns (dx/dt, jacobian) at a given state.
/// Uses RK4 integration with periodic QR reorthonormalization.
pub fn lyapunov_exponents_flow<F>(
    f: F,
    x0: &DVector<f64>,
    dt: f64,
    n_steps: usize,
    renorm_interval: usize,
) -> Vec<f64>
where
    F: Fn(&DVector<f64>) -> (DVector<f64>, DMatrix<f64>),
{
    let dim = x0.nrows();
    let mut x = x0.clone();
    let mut q = DMatrix::<f64>::identity(dim, dim);
    let mut sums = vec![0.0f64; dim];
    let mut renorm_count = 0usize;

    // Transient
    let transient = n_steps / 10;
    for _ in 0..transient {
        let (dx, _) = f(&x);
        x += &(dx * dt);
    }

    for step in 0..n_steps {
        // Integrate the variational equations using RK4
        let (dx1, jac1) = f(&x);
        let x2 = &x + &(dx1.clone() * (dt * 0.5));
        let (dx2, jac2) = f(&x2);
        let x3 = &x + &(dx2.clone() * (dt * 0.5));
        let (dx3, jac3) = f(&x3);
        let x4 = &x + &(dx3.clone() * dt);
        let (dx4, jac4) = f(&x4);

        // Update state with RK4
        let dx = &(dx1 * 1.0) + &(dx2 * 2.0) + &(dx3 * 2.0) + &dx4;
        x += &(dx * (dt / 6.0));

        // Average Jacobian for variational update
        let avg_jac = &jac1 + &jac2 * 2.0 + &jac3 * 2.0 + &jac4;
        let avg_jac = avg_jac * (1.0 / 6.0);

        // Update tangent vectors: dq/dt = J * q
        let dq = &avg_jac * &q;
        q += &(dq * dt);

        // Reorthonormalize
        if (step + 1) % renorm_interval == 0 {
            let qr = q.qr();
            let r = qr.r();
            for j in 0..dim {
                let diag = r[(j, j)].abs();
                if diag > 0.0 {
                    sums[j] += diag.ln();
                }
            }
            q = qr.q();
            renorm_count += 1;
        }
    }

    let total_time = renorm_count as f64 * (renorm_interval as f64 * dt);
    sums.iter().map(|s| s / total_time).collect()
}

/// Maximum Lyapunov exponent (convenience wrapper).
pub fn max_lyapunov(exponents: &[f64]) -> f64 {
    exponents.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

/// Check if a system exhibits chaos (positive max Lyapunov exponent).
pub fn is_chaotic(exponents: &[f64]) -> bool {
    max_lyapunov(exponents) > 0.0
}

/// Kaplan-Yorke dimension estimate from Lyapunov spectrum.
pub fn kaplan_yorke_dimension(exponents: &[f64]) -> f64 {
    let mut sorted: Vec<f64> = exponents.to_vec();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let mut sum = 0.0f64;
    let mut j = 0usize;
    for (i, &exp) in sorted.iter().enumerate() {
        let new_sum = sum + exp;
        if new_sum >= 0.0 {
            sum = new_sum;
            j = i;
        } else {
            break;
        }
    }
    if j + 1 < sorted.len() {
        let d_ky = (j as f64) + 1.0 + sum / sorted[j + 1].abs();
        d_ky.max(0.0)
    } else {
        j as f64 + 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DVector, DMatrix};

    #[test]
    fn test_logistic_lyapunov_stable() {
        let f = |x: f64| 3.2 * x * (1.0 - x);
        let df = |x: f64| 3.2 * (1.0 - 2.0 * x);
        let le = lyapunov_exponent_map(f, df, 0.1, 10000);
        assert!(le < 0.0, "Period-2 should have negative LE, got {}", le);
    }

    #[test]
    fn test_logistic_lyapunov_chaotic() {
        let f = |x: f64| 4.0 * x * (1.0 - x);
        let df = |x: f64| 4.0 * (1.0 - 2.0 * x);
        let le = lyapunov_exponent_map(f, df, 0.1, 50000);
        assert!((le - 2.0f64.ln()).abs() < 0.1,
            "LE should be ~ln(2)={:.3}, got {:.3}", 2.0f64.ln(), le);
    }

    #[test]
    fn test_logistic_lyapunov_r3() {
        let f = |x: f64| 3.0 * x * (1.0 - x);
        let df = |x: f64| 3.0 * (1.0 - 2.0 * x);
        let le = lyapunov_exponent_map(f, df, 0.2, 50000);
        assert!(le.abs() < 0.05, "LE near 0 at r=3, got {}", le);
    }

    #[test]
    fn test_tent_map_lyapunov() {
        let f = |x: f64| if x < 0.5 { 2.0 * x } else { 2.0 * (1.0 - x) };
        let df = |x: f64| if x < 0.5 { 2.0 } else { -2.0 };
        let le = lyapunov_exponent_map(f, df, 0.3, 50000);
        assert!((le - 2.0f64.ln()).abs() < 0.1,
            "Tent map LE should be ~ln(2), got {:.3}", le);
    }

    #[test]
    fn test_henon_lyapunov_correct() {
        let f = |x: &DVector<f64>| -> (DVector<f64>, DMatrix<f64>) {
            let xn = 1.0 - 1.4 * x[0] * x[0] + 0.3 * x[1];
            let yn = x[0];
            let new_x = DVector::from_vec(vec![xn, yn]);
            let jac = DMatrix::from_row_slice(2, 2, &[
                -2.8 * x[0], 0.3,
                1.0, 0.0,
            ]);
            (new_x, jac)
        };
        let g = |x: &DVector<f64>| -> (DVector<f64>, DMatrix<f64>) {
            f(x)
        };
        let x0 = DVector::from_vec(vec![0.0, 0.0]);
        let exps = lyapunov_exponents_map(f, g, &x0, 20000);
        assert!(exps[0] > 0.0, "Hénon max LE should be positive, got {}", exps[0]);
        assert!(exps[1] < 0.0, "Hénon min LE should be negative, got {}", exps[1]);
    }

    #[test]
    fn test_lorenz_lyapunov() {
        let sigma = 10.0_f64;
        let rho = 28.0_f64;
        let beta = 8.0 / 3.0;
        let f = move |x: &DVector<f64>| -> (DVector<f64>, DMatrix<f64>) {
            let dx = DVector::from_vec(vec![
                sigma * (x[1] - x[0]),
                x[0] * (rho - x[2]) - x[1],
                x[0] * x[1] - beta * x[2],
            ]);
            let jac = DMatrix::from_row_slice(3, 3, &[
                -sigma, sigma, 0.0,
                rho - x[2], -1.0, -x[0],
                x[1], x[0], -beta,
            ]);
            (dx, jac)
        };
        let x0 = DVector::from_vec(vec![1.0, 1.0, 1.0]);
        let exps = lyapunov_exponents_flow(f, &x0, 0.01, 50000, 10);
        assert!(exps[0] > 0.3, "Lorenz max LE should be >0.3, got {}", exps[0]);
        assert!(is_chaotic(&exps), "Lorenz should be chaotic");
        let dky = kaplan_yorke_dimension(&exps);
        assert!(dky > 1.8 && dky < 2.5, "KY dimension should be ~2.06, got {}", dky);
    }

    #[test]
    fn test_kaplan_yorke_dimension() {
        let exps = vec![0.906, 0.0, -14.572];
        let dky = kaplan_yorke_dimension(&exps);
        assert!((dky - 2.062).abs() < 0.01, "KY dim should be ~2.062, got {}", dky);
    }

    #[test]
    fn test_is_chaotic() {
        assert!(is_chaotic(&[0.5, -2.0]));
        assert!(!is_chaotic(&[-0.5, -2.0]));
    }
}
