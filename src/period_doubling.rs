//! Period-doubling cascade and Feigenbaum constants.

/// First Feigenbaum constant: δ ≈ 4.6692...
/// Ratio between successive bifurcation parameter intervals.
pub const FEIGENBAUM_DELTA: f64 = 4.66920160910299067185;

/// Second Feigenbaum constant: α ≈ 2.5029...
/// Ratio between successive period-doubling scales in state space.
pub const FEIGENBAUM_ALPHA: f64 = 2.50290787509589282228;

/// Find the period-doubling bifurcation points of the logistic map.
///
/// Returns the parameter values where period doubling occurs.
pub fn logistic_bifurcation_points(r_min: f64, r_max: f64, resolution: usize) -> Vec<f64> {
    let dr = (r_max - r_min) / resolution as f64;
    let mut points = Vec::new();
    let mut prev_period = 0usize;

    for i in 0..=resolution {
        let r = r_min + dr * i as f64;
        let period = detect_period(r, 1000, 200);
        if prev_period > 0 && period != prev_period && period > prev_period {
            points.push(r - dr * 0.5);
        }
        prev_period = period;
    }

    points
}

/// Detect the period of the logistic map orbit at parameter r.
fn detect_period(r: f64, transient: usize, test: usize) -> usize {
    let mut x = 0.5_f64;
    for _ in 0..transient {
        x = r * x * (1.0 - x);
        if !x.is_finite() { return 0; }
    }

    let x0 = x;
    for n in 1..=test {
        x = r * x * (1.0 - x);
        if !x.is_finite() { return 0; }
        if (x - x0).abs() < 1e-8 {
            return n;
        }
    }
    0 // Could not determine period (chaotic or very long period)
}

/// Estimate the Feigenbaum δ from a series of bifurcation points.
///
/// δ_n = (r_n - r_{n-1}) / (r_{n+1} - r_n)
pub fn estimate_feigenbaum_delta(bif_points: &[f64]) -> f64 {
    if bif_points.len() < 3 { return 0.0; }
    let n = bif_points.len() - 2;
    let rn = bif_points[n];
    let rn_1 = bif_points[n - 1];
    let rn1 = bif_points[n + 1];
    (rn_1 - rn) / (rn1 - rn)
}

/// Verify Feigenbaum universality: compute δ for the logistic map at known bifurcation points.
pub fn verify_feigenbaum_delta() -> f64 {
    // Known logistic map bifurcation points
    let _r1 = 3.0;        // period 1 -> 2
    let r2 = 3.44949;    // period 2 -> 4
    let r3 = 3.54409;    // period 4 -> 8
    let r4 = 3.564407;   // period 8 -> 16

    // δ from r2, r3, r4
    let delta = (r3 - r2) / (r4 - r3);
    delta
}

/// Compute the superstable orbit parameter for period 2^n in the logistic map.
///
/// A superstable orbit occurs when the critical point x=0.5 maps back to itself
/// after 2^n iterations.
pub fn superstable_parameter(period_power: u32) -> f64 {
    let target_period = 1usize << period_power;
    let r_min = 2.5;
    let r_max = 4.0;
    let steps = 100000;
    let dr = (r_max - r_min) / steps as f64;

    let mut best_r = 0.0;
    let mut best_err = f64::INFINITY;

    for i in 0..=steps {
        let r = r_min + dr * i as f64;
        let err = superstable_error(r, target_period);
        if err < best_err {
            best_err = err;
            best_r = r;
        }
    }

    best_r
}

fn superstable_error(r: f64, period: usize) -> f64 {
    let mut x = 0.5;
    for _ in 0..period {
        x = r * x * (1.0 - x);
    }
    (x - 0.5).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feigenbaum_delta_constant() {
        assert!((FEIGENBAUM_DELTA - 4.669).abs() < 0.001);
    }

    #[test]
    fn test_feigenbaum_alpha_constant() {
        assert!((FEIGENBAUM_ALPHA - 2.503).abs() < 0.001);
    }

    #[test]
    fn test_verify_feigenbaum_delta() {
        let delta = verify_feigenbaum_delta();
        // Should be approaching 4.669...
        assert!((delta - FEIGENBAUM_DELTA).abs() < 0.5,
            "δ should approach 4.669, got {:.3}", delta);
        assert!(delta > 3.0, "δ should be > 3, got {:.3}", delta);
    }

    #[test]
    fn test_logistic_bifurcation_detection() {
        let points = logistic_bifurcation_points(2.5, 3.6, 10000);
        assert!(!points.is_empty(), "Should find bifurcation points");
        // First bifurcation should be near r=3
        assert!(points[0] > 2.9 && points[0] < 3.1,
            "First bifurcation near r=3, got {:.3}", points[0]);
    }

    #[test]
    fn test_period_detection() {
        // r=2.5: fixed point (period 1)
        assert_eq!(detect_period(2.5, 1000, 200), 1);
        // r=3.2: period 2
        assert_eq!(detect_period(3.2, 1000, 200), 2);
        // r=3.5: period 4
        assert_eq!(detect_period(3.5, 1000, 200), 4);
    }

    #[test]
    fn test_estimate_feigenbaum() {
        let points = vec![3.0, 3.44949, 3.54409, 3.564407];
        let delta = estimate_feigenbaum_delta(&points);
        assert!(delta.abs() > 3.0 && delta.abs() < 7.0, "δ estimate should be reasonable, got {:.3}", delta);
    }
}
