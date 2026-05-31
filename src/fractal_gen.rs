//! Fractal generation: Mandelbrot set, Julia sets, escape time algorithm.

use serde::{Serialize, Deserialize};

/// Result of a fractal membership test.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FractalResult {
    /// Whether the point is in the set.
    pub in_set: bool,
    /// Number of iterations before escape (None if didn't escape).
    pub iterations: Option<usize>,
    /// Smooth escape time for coloring.
    pub smooth_escape: Option<f64>,
}

/// Test if a point c is in the Mandelbrot set.
///
/// Iterates z_{n+1} = z_n^2 + c starting from z_0 = 0.
/// Returns after `max_iter` iterations if the orbit hasn't escaped.
pub fn mandelbrot(c_re: f64, c_im: f64, max_iter: usize) -> FractalResult {
    let mut z_re = 0.0_f64;
    let mut z_im = 0.0_f64;

    for i in 0..max_iter {
        let z_re2 = z_re * z_re;
        let z_im2 = z_im * z_im;
        if z_re2 + z_im2 > 4.0 {
            // Escaped — compute smooth iteration count
            let log_zn = (z_re2 + z_im2).ln() / 2.0;
            let nu = (log_zn / 2.0_f64.ln()).ln();
            let smooth = i as f64 + 1.0 - nu;
            return FractalResult {
                in_set: false,
                iterations: Some(i),
                smooth_escape: Some(smooth),
            };
        }
        z_im = 2.0 * z_re * z_im + c_im;
        z_re = z_re2 - z_im2 + c_re;
    }

    FractalResult {
        in_set: true,
        iterations: None,
        smooth_escape: None,
    }
}

/// Test if a point z_0 is in the Julia set for parameter c.
///
/// Iterates z_{n+1} = z_n^2 + c.
pub fn julia(z_re: f64, z_im: f64, c_re: f64, c_im: f64, max_iter: usize) -> FractalResult {
    let mut zr = z_re;
    let mut zi = z_im;

    for i in 0..max_iter {
        let zr2 = zr * zr;
        let zi2 = zi * zi;
        if zr2 + zi2 > 4.0 {
            let log_zn = (zr2 + zi2).ln() / 2.0;
            let nu = (log_zn / 2.0_f64.ln()).ln();
            let smooth = i as f64 + 1.0 - nu;
            return FractalResult {
                in_set: false,
                iterations: Some(i),
                smooth_escape: Some(smooth),
            };
        }
        let new_zr = zr2 - zi2 + c_re;
        zi = 2.0 * zr * zi + c_im;
        zr = new_zr;
    }

    FractalResult {
        in_set: true,
        iterations: None,
        smooth_escape: None,
    }
}

/// Compute a grid of Mandelbrot escape times.
pub fn mandelbrot_grid(
    x_range: (f64, f64),
    y_range: (f64, f64),
    width: usize,
    height: usize,
    max_iter: usize,
) -> Vec<Vec<Option<usize>>> {
    let mut grid = Vec::with_capacity(height);
    for row in 0..height {
        let y = y_range.0 + (y_range.1 - y_range.0) * row as f64 / (height - 1) as f64;
        let mut row_data = Vec::with_capacity(width);
        for col in 0..width {
            let x = x_range.0 + (x_range.1 - x_range.0) * col as f64 / (width - 1) as f64;
            let result = mandelbrot(x, y, max_iter);
            row_data.push(result.iterations);
        }
        grid.push(row_data);
    }
    grid
}

/// Compute a grid of Julia set escape times.
pub fn julia_grid(
    c_re: f64,
    c_im: f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    width: usize,
    height: usize,
    max_iter: usize,
) -> Vec<Vec<Option<usize>>> {
    let mut grid = Vec::with_capacity(height);
    for row in 0..height {
        let y = y_range.0 + (y_range.1 - y_range.0) * row as f64 / (height - 1) as f64;
        let mut row_data = Vec::with_capacity(width);
        for col in 0..width {
            let x = x_range.0 + (x_range.1 - x_range.0) * col as f64 / (width - 1) as f64;
            let result = julia(x, y, c_re, c_im, max_iter);
            row_data.push(result.iterations);
        }
        grid.push(row_data);
    }
    grid
}

/// Burning Ship fractal: z_{n+1} = (|Re(z_n)| + i|Im(z_n)|)^2 + c
pub fn burning_ship(c_re: f64, c_im: f64, max_iter: usize) -> FractalResult {
    let mut z_re = 0.0_f64;
    let mut z_im = 0.0_f64;

    for i in 0..max_iter {
        let z_re2 = z_re * z_re;
        let z_im2 = z_im * z_im;
        if z_re2 + z_im2 > 4.0 {
            let log_zn = (z_re2 + z_im2).ln() / 2.0;
            let nu = (log_zn / 2.0_f64.ln()).ln();
            let smooth = i as f64 + 1.0 - nu;
            return FractalResult {
                in_set: false,
                iterations: Some(i),
                smooth_escape: Some(smooth),
            };
        }
        z_im = (2.0 * z_re * z_im).abs() + c_im;
        z_re = z_re2 - z_im2 + c_re;
    }

    FractalResult {
        in_set: true,
        iterations: None,
        smooth_escape: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mandelbrot_center_in_set() {
        // c = 0 is clearly in the Mandelbrot set
        let r = mandelbrot(0.0, 0.0, 100);
        assert!(r.in_set);
    }

    #[test]
    fn test_mandelbrot_escape() {
        // c = 2 + 0i escapes immediately
        let r = mandelbrot(2.0, 0.0, 100);
        assert!(!r.in_set);
        assert!(r.iterations.unwrap() <= 2);
    }

    #[test]
    fn test_mandelbrot_main_cardioid() {
        // Points inside the main cardioid
        let r = mandelbrot(-0.5, 0.0, 100);
        assert!(r.in_set);
        let r2 = mandelbrot(0.25, 0.0, 100);
        assert!(r2.in_set);
    }

    #[test]
    fn test_mandelbrot_period2_bulb() {
        // c = -1 is in the period-2 bulb
        let r = mandelbrot(-1.0, 0.0, 200);
        assert!(r.in_set);
    }

    #[test]
    fn test_mandelbrot_outside() {
        // c = 0.5 + 0.5i is outside
        let r = mandelbrot(0.5, 0.5, 100);
        assert!(!r.in_set);
    }

    #[test]
    fn test_julia_in_set() {
        // c = 0: z -> z^2, points with |z| <= 1 stay bounded
        let r = julia(0.5, 0.0, 0.0, 0.0, 100);
        assert!(r.in_set);
    }

    #[test]
    fn test_julia_escape() {
        let r = julia(2.0, 0.0, 0.0, 0.0, 100);
        assert!(!r.in_set);
    }

    #[test]
    fn test_julia_dendrite() {
        // c = i gives a dendrite Julia set; z=0 should be in set
        let r = julia(0.0, 0.0, 0.0, 1.0, 200);
        assert!(r.in_set);
    }

    #[test]
    fn test_julia_cantor_dust() {
        // c = -2 gives a line segment Julia set; z=0 is in set
        let r = julia(0.0, 0.0, -2.0, 0.0, 100);
        assert!(r.in_set);
    }

    #[test]
    fn test_mandelbrot_grid() {
        let grid = mandelbrot_grid((-2.0, 1.0), (-1.5, 1.5), 100, 100, 50);
        assert_eq!(grid.len(), 100);
        assert_eq!(grid[0].len(), 100);
        // Center point should be in set
        assert!(grid[50][35].is_none()); // in set = None iterations
    }

    #[test]
    fn test_burning_ship() {
        // c = 0 should be in set
        let r = burning_ship(0.0, 0.0, 100);
        assert!(r.in_set);
        // c = 2 should escape
        let r = burning_ship(2.0, 0.0, 100);
        assert!(!r.in_set);
    }
}
