//! Agent unpredictability analysis — measuring and characterizing chaotic behavior in agent systems.

use serde::{Serialize, Deserialize};
// Sensitivity and fractal dimension imports used by analysis functions
use crate::symbolic::{Partition, encode_itinerary, topological_entropy, shannon_entropy_rate};

/// Metric type for agent behavior.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BehaviorMetric {
    /// Raw scalar observations.
    Scalar(Vec<f64>),
    /// Vector observations (e.g., multi-dimensional state).
    Vector(Vec<Vec<f64>>),
    /// Categorical/discrete observations.
    Categorical(Vec<usize>),
}

/// Result of unpredictability analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnpredictabilityReport {
    /// Estimated Lyapunov exponent (positive = chaotic).
    pub lyapunov_exponent: f64,
    /// Whether the behavior is classified as chaotic.
    pub is_chaotic: bool,
    /// Correlation dimension estimate (if applicable).
    pub correlation_dim: Option<f64>,
    /// Shannon entropy rate of symbolic encoding.
    pub entropy_rate: f64,
    /// Topological entropy estimate.
    pub topological_entropy: f64,
    /// Predictability score (0 = perfectly predictable, 1 = maximally unpredictable).
    pub unpredictability_score: f64,
    /// Qualitative description.
    pub description: String,
}

/// Analyze unpredictability of a scalar time series.
///
/// Fits a local 1D map model and computes Lyapunov exponent,
/// then supplements with entropy-based measures.
pub fn analyze_scalar_unpredictability(
    series: &[f64],
    partition_boundaries: Vec<f64>,
) -> UnpredictabilityReport {
    if series.len() < 10 {
        return UnpredictabilityReport {
            lyapunov_exponent: 0.0,
            is_chaotic: false,
            correlation_dim: None,
            entropy_rate: 0.0,
            topological_entropy: 0.0,
            unpredictability_score: 0.0,
            description: "Insufficient data".to_string(),
        };
    }

    // Estimate Lyapunov exponent using nearest-neighbor divergence
    let le = estimate_lyapunov_from_series(series);

    // Symbolic encoding
    let partition = Partition::new(partition_boundaries);
    let symbols = encode_itinerary(series, &partition);

    // Entropy measures
    let entropy_rate = shannon_entropy_rate(&symbols);
    let topo_entropy = topological_entropy(&symbols, 5.min(series.len() / 2));

    // Predictability score: combination of LE and entropy
    let le_component = (le / 2.0).min(1.0).max(0.0);
    let entropy_component = (entropy_rate / 2.0f64.ln()).min(1.0).max(0.0);
    let unpredictability = 0.5 * le_component + 0.5 * entropy_component;

    let chaotic = le > 0.0;
    let description = if chaotic {
        "Chaotic behavior detected — sensitive dependence on initial conditions".to_string()
    } else if le > -0.5 {
        "Borderline behavior — near transition to chaos".to_string()
    } else {
        "Regular/predictable behavior".to_string()
    };

    UnpredictabilityReport {
        lyapunov_exponent: le,
        is_chaotic: chaotic,
        correlation_dim: None,
        entropy_rate: entropy_rate,
        topological_entropy: topo_entropy,
        unpredictability_score: unpredictability,
        description,
    }
}

/// Estimate the maximum Lyapunov exponent from a scalar time series
/// using the Rosenstein et al. method (nearest-neighbor divergence).
fn estimate_lyapunov_from_series(series: &[f64]) -> f64 {
    let n = series.len();
    if n < 20 { return 0.0; }

    let mut min_dist = f64::INFINITY;
    for i in 1..n {
        let d = (series[i] - series[i - 1]).abs();
        if d > 0.0 && d < min_dist {
            min_dist = d;
        }
    }

    // Estimate local divergence rate
    let mut divergences = Vec::new();
    let max_lag = n / 4;

    for i in 0..n.saturating_sub(max_lag) {
        // Find nearest neighbor
        let mut best_j = None;
        let mut best_dist = f64::INFINITY;
        for j in 0..n.saturating_sub(max_lag) {
            if (j as isize - i as isize).unsigned_abs() < 5 { continue; }
            let d = (series[i] - series[j]).abs();
            if d < best_dist && d > 0.0 {
                best_dist = d;
                best_j = Some(j);
            }
        }

        if let Some(j) = best_j {
            let mut sum = 0.0;
            let mut count = 0;
            for k in 1..max_lag {
                if i + k < n && j + k < n {
                    let dist = (series[i + k] - series[j + k]).abs();
                    if dist > 0.0 {
                        sum += dist.ln();
                        count += 1;
                    }
                }
            }
            if count > 0 {
                divergences.push(sum / count as f64);
            }
        }
    }

    if divergences.is_empty() {
        return 0.0;
    }

    let avg_div: f64 = divergences.iter().sum::<f64>() / divergences.len() as f64;
    avg_div / max_lag as f64
}

/// Agent behavior profile — summary statistics for monitoring.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentProfile {
    pub variance: f64,
    pub autocorrelation_lag1: f64,
    pub entropy_rate: f64,
    pub unpredictability_score: f64,
    pub behavioral_regime: BehavioralRegime,
}

/// Classification of agent behavioral regime.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum BehavioralRegime {
    FixedPoint,
    Periodic,
    Quasiperiodic,
    Chaotic,
    Unknown,
}

/// Build an agent behavior profile from a scalar time series.
pub fn build_agent_profile(series: &[f64]) -> AgentProfile {
    let n = series.len() as f64;
    if n < 3.0 {
        return AgentProfile {
            variance: 0.0,
            autocorrelation_lag1: 0.0,
            entropy_rate: 0.0,
            unpredictability_score: 0.0,
            behavioral_regime: BehavioralRegime::Unknown,
        };
    }

    // Variance
    let mean: f64 = series.iter().sum::<f64>() / n;
    let variance: f64 = series.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;

    // Autocorrelation at lag 1
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 1..series.len() {
        num += (series[i] - mean) * (series[i - 1] - mean);
    }
    for x in series {
        den += (x - mean).powi(2);
    }
    let acf1 = if den.abs() > 1e-15 { num / den } else { 0.0 };

    // Determine regime
    let min_val = series.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_val = series.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max_val - min_val;

    let report = analyze_scalar_unpredictability(series, vec![mean]);

    let regime = if variance < 1e-10 {
        BehavioralRegime::FixedPoint
    } else if !report.is_chaotic && acf1.abs() > 0.8 {
        BehavioralRegime::Periodic
    } else if report.is_chaotic {
        BehavioralRegime::Chaotic
    } else if acf1.abs() < 0.3 && range > 0.0 {
        BehavioralRegime::Quasiperiodic
    } else {
        BehavioralRegime::Unknown
    };

    AgentProfile {
        variance,
        autocorrelation_lag1: acf1,
        entropy_rate: report.entropy_rate,
        unpredictability_score: report.unpredictability_score,
        behavioral_regime: regime,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_series() {
        let series = vec![1.0; 100];
        let profile = build_agent_profile(&series);
        assert_eq!(profile.behavioral_regime, BehavioralRegime::FixedPoint);
        assert!(profile.variance < 1e-10);
    }

    #[test]
    fn test_periodic_series() {
        let series: Vec<f64> = (0..200).map(|i| (i as f64 * 0.1).sin()).collect();
        let profile = build_agent_profile(&series);
        assert!(profile.autocorrelation_lag1.abs() > 0.5);
    }

    #[test]
    fn test_chaotic_logistic_series() {
        let mut x = 0.1_f64;
        let mut series = Vec::new();
        for _ in 0..500 {
            series.push(x);
            x = 4.0 * x * (1.0 - x);
        }
        let report = analyze_scalar_unpredictability(&series, vec![0.5]);
        assert!(report.lyapunov_exponent > -0.5, "Should show some divergence, got LE={}", report.lyapunov_exponent);
    }

    #[test]
    fn test_unpredictability_report() {
        let series: Vec<f64> = (0..100).map(|i| (i as f64 * 0.5).sin()).collect();
        let report = analyze_scalar_unpredictability(&series, vec![0.0]);
        assert!(!report.is_chaotic);
        assert!(!report.description.is_empty());
    }

    #[test]
    fn test_insufficient_data() {
        let series = vec![1.0, 2.0];
        let report = analyze_scalar_unpredictability(&series, vec![1.5]);
        assert_eq!(report.unpredictability_score, 0.0);
    }

    #[test]
    fn test_agent_profile_random() {
        // Pseudo-random series
        let series: Vec<f64> = (0..500).map(|i| {
            let x = ((i as f64 * 127.1 + 311.7).sin() * 43758.5453) % 1.0;
            x.abs()
        }).collect();
        let profile = build_agent_profile(&series);
        assert!(profile.variance > 0.01, "Random series should have variance > 0");
    }
}
