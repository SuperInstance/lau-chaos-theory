//! Fractal dimension estimation: box-counting, correlation, information dimension.

use std::collections::HashMap;

/// Box-counting dimension estimate.
///
/// Given a set of points, estimates D_box = lim( ln N(ε) / ln(1/ε) )
/// Uses multiple box sizes and linear regression on log-log scale.
pub fn box_counting_dimension(points: &[(f64, f64)], scales: &[f64]) -> f64 {
    let mut log_inv_scale: Vec<f64> = Vec::new();
    let mut log_count: Vec<f64> = Vec::new();

    for &eps in scales {
        if eps <= 0.0 { continue; }
        let mut boxes: HashMap<(i64, i64), ()> = HashMap::new();
        for &(x, y) in points {
            let bx = (x / eps).floor() as i64;
            let by = (y / eps).floor() as i64;
            boxes.insert((bx, by), ());
        }
        let n = boxes.len() as f64;
        if n > 0.0 {
            log_inv_scale.push((1.0 / eps).ln());
            log_count.push(n.ln());
        }
    }

    if log_inv_scale.len() < 2 {
        return 0.0;
    }

    linear_regression(&log_inv_scale, &log_count)
}

/// Correlation dimension (Grassberger-Procaccia algorithm).
///
/// C(ε) = (2 / N(N-1)) Σ_{i<j} Θ(ε - |x_i - x_j|)
/// D_corr = lim ln C(ε) / ln(ε)
pub fn correlation_dimension(points: &[(f64, f64)], scales: &[f64]) -> f64 {
    let n = points.len();
    if n < 2 { return 0.0; }

    let mut log_eps: Vec<f64> = Vec::new();
    let mut log_corr: Vec<f64> = Vec::new();

    let total_pairs = (n * (n - 1)) as f64 / 2.0;

    for &eps in scales {
        if eps <= 0.0 { continue; }
        let mut count = 0usize;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = points[i].0 - points[j].0;
                let dy = points[i].1 - points[j].1;
                if (dx * dx + dy * dy).sqrt() < eps {
                    count += 1;
                }
            }
        }
        let c = count as f64 / total_pairs;
        if c > 0.0 {
            log_eps.push(eps.ln());
            log_corr.push(c.ln());
        }
    }

    if log_eps.len() < 2 {
        return 0.0;
    }

    linear_regression(&log_eps, &log_corr)
}

/// Information dimension.
///
/// D_1 = lim Σ p_i ln(1/p_i) / ln(1/ε)
/// where p_i is the probability of a point landing in box i.
pub fn information_dimension(points: &[(f64, f64)], scales: &[f64]) -> f64 {
    let n = points.len() as f64;

    let mut log_inv_scale: Vec<f64> = Vec::new();
    let mut entropy: Vec<f64> = Vec::new();

    for &eps in scales {
        if eps <= 0.0 { continue; }
        let mut box_count: HashMap<(i64, i64), usize> = HashMap::new();
        for &(x, y) in points {
            let bx = (x / eps).floor() as i64;
            let by = (y / eps).floor() as i64;
            *box_count.entry((bx, by)).or_insert(0) += 1;
        }
        let h: f64 = box_count.values()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / n;
                -p * p.ln()
            })
            .sum();

        if h > 0.0 {
            log_inv_scale.push((1.0 / eps).ln());
            entropy.push(h);
        }
    }

    if log_inv_scale.len() < 2 {
        return 0.0;
    }

    linear_regression(&log_inv_scale, &entropy)
}

/// Simple linear regression: returns slope of y ~ x.
fn linear_regression(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx: f64 = x.iter().sum::<f64>() / n;
    let my: f64 = y.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        num += dx * (y[i] - my);
        den += dx * dx;
    }
    if den.abs() < 1e-15 { 0.0 } else { num / den }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate points on a line segment — should give dimension ~1.
    fn line_points(n: usize) -> Vec<(f64, f64)> {
        (0..n).map(|i| {
            let t = i as f64 / n as f64;
            (t, 0.5 * t)
        }).collect()
    }

    /// Generate points uniformly filling a square — should give dimension ~2.
    fn square_points(n: usize) -> Vec<(f64, f64)> {
        let side = (n as f64).sqrt() as usize;
        (0..side * side).map(|i| {
            let row = i / side;
            let col = i % side;
            (row as f64 / side as f64 + 0.5 / side as f64, 
             col as f64 / side as f64 + 0.5 / side as f64)
        }).collect()
    }

    #[test]
    fn test_box_counting_line() {
        let pts = line_points(500);
        let scales: Vec<f64> = (1..8).map(|i| 0.5_f64.powi(i)).collect();
        let dim = box_counting_dimension(&pts, &scales);
        assert!((dim - 1.0).abs() < 0.15, "Line should have dim ~1, got {}", dim);
    }

    #[test]
    fn test_box_counting_square() {
        let pts = square_points(2500);
        let scales: Vec<f64> = (1..6).map(|i| 0.5_f64.powi(i)).collect();
        let dim = box_counting_dimension(&pts, &scales);
        assert!((dim - 2.0).abs() < 0.3, "Square should have dim ~2, got {}", dim);
    }

    #[test]
    fn test_correlation_dimension_line() {
        let pts = line_points(200);
        let scales: Vec<f64> = vec![0.5, 0.3, 0.2, 0.1, 0.05, 0.02];
        let dim = correlation_dimension(&pts, &scales);
        assert!((dim - 1.0).abs() < 0.3, "Line correlation dim ~1, got {}", dim);
    }

    #[test]
    fn test_correlation_dimension_square() {
        let pts = square_points(400);
        let scales: Vec<f64> = vec![0.5, 0.3, 0.2, 0.1, 0.05, 0.02];
        let dim = correlation_dimension(&pts, &scales);
        assert!((dim - 2.0).abs() < 0.3, "Square correlation dim ~2, got {}", dim);
    }

    #[test]
    fn test_information_dimension_line() {
        let pts = line_points(300);
        let scales: Vec<f64> = vec![0.5, 0.3, 0.2, 0.1, 0.05];
        let dim = information_dimension(&pts, &scales);
        assert!((dim - 1.0).abs() < 0.2, "Line info dim ~1, got {}", dim);
    }

    #[test]
    fn test_information_dimension_square() {
        let pts = square_points(2500);
        let scales: Vec<f64> = vec![0.5, 0.3, 0.2, 0.1, 0.05];
        let dim = information_dimension(&pts, &scales);
        assert!((dim - 2.0).abs() < 0.4, "Square info dim ~2, got {}", dim);
    }
}
