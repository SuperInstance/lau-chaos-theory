//! Symbolic dynamics: shift spaces, itinerary encoding, topological entropy.

use serde::{Serialize, Deserialize};

/// A partition of the state space into labeled regions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Partition {
    /// Boundaries that divide the state space (sorted).
    pub boundaries: Vec<f64>,
}

impl Partition {
    /// Create a new partition from sorted boundaries.
    pub fn new(boundaries: Vec<f64>) -> Self {
        let mut b = boundaries;
        b.sort_by(|a, b| a.partial_cmp(b).unwrap());
        Self { boundaries: b }
    }

    /// Get the symbol (region index) for a given value.
    pub fn symbol(&self, x: f64) -> usize {
        let mut idx = 0;
        for &b in &self.boundaries {
            if x >= b {
                idx += 1;
            } else {
                break;
            }
        }
        idx
    }

    /// Number of symbols.
    pub fn n_symbols(&self) -> usize {
        self.boundaries.len() + 1
    }
}

/// Encode a trajectory into a symbolic sequence using a given partition.
pub fn encode_itinerary(trajectory: &[f64], partition: &Partition) -> Vec<usize> {
    trajectory.iter().map(|&x| partition.symbol(x)).collect()
}

/// Full shift space of n symbols.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShiftSpace {
    pub n_symbols: usize,
    /// Forbidden transitions: (from_symbol, to_symbol)
    pub forbidden: Vec<(usize, usize)>,
}

impl ShiftSpace {
    /// Create a full shift on n symbols (no forbidden transitions).
    pub fn full(n: usize) -> Self {
        Self { n_symbols: n, forbidden: vec![] }
    }

    /// Check if a transition is allowed.
    pub fn is_allowed(&self, from: usize, to: usize) -> bool {
        !self.forbidden.contains(&(from, to))
    }

    /// Count allowed transitions for a given symbol.
    pub fn allowed_transitions(&self, from: usize) -> usize {
        (0..self.n_symbols)
            .filter(|&to| self.is_allowed(from, to))
            .count()
    }
}

/// Estimate topological entropy from a symbolic sequence.
///
/// Uses the formula: h_top ≈ lim ln(N(n)) / n where N(n) is the number of
/// distinct n-blocks. We estimate from finite data.
pub fn topological_entropy(sequence: &[usize], block_size: usize) -> f64 {
    if sequence.len() < block_size { return 0.0; }

    let mut blocks = std::collections::HashSet::new();
    for i in 0..=(sequence.len() - block_size) {
        let block: Vec<usize> = sequence[i..(i + block_size)].to_vec();
        blocks.insert(block);
    }

    let n_blocks = blocks.len() as f64;
    if n_blocks <= 1.0 { return 0.0; }
    n_blocks.ln() / block_size as f64
}

/// Compute the Shannon entropy rate from a symbolic sequence.
pub fn shannon_entropy_rate(sequence: &[usize]) -> f64 {
    use std::collections::HashMap;
    if sequence.len() < 2 { return 0.0; }

    let mut transitions: HashMap<(usize, usize), usize> = HashMap::new();
    let mut from_counts: HashMap<usize, usize> = HashMap::new();

    for i in 0..(sequence.len() - 1) {
        let from = sequence[i];
        let to = sequence[i + 1];
        *transitions.entry((from, to)).or_insert(0) += 1;
        *from_counts.entry(from).or_insert(0) += 1;
    }

    let mut h = 0.0;
    for ((from, _to), count) in &transitions {
        let total = from_counts[from];
        if *count > 0 && total > 0 {
            let p = *count as f64 / total as f64;
            h -= p * p.ln();
        }
    }

    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_symbol() {
        let p = Partition::new(vec![0.0, 1.0, 2.0]);
        assert_eq!(p.symbol(-1.0), 0);
        assert_eq!(p.symbol(0.0), 1);
        assert_eq!(p.symbol(0.5), 1);
        assert_eq!(p.symbol(1.0), 2);
        assert_eq!(p.symbol(1.5), 2);
        assert_eq!(p.symbol(2.0), 3);
        assert_eq!(p.symbol(3.0), 3);
        assert_eq!(p.n_symbols(), 4);
    }

    #[test]
    fn test_encode_logistic_map() {
        // Logistic map r=4, partition at x=0.5
        let p = Partition::new(vec![0.5]);
        let f = |x: f64| 4.0 * x * (1.0 - x);
        let mut x = 0.2_f64;
        let mut traj = Vec::new();
        for _ in 0..100 {
            traj.push(x);
            x = f(x);
        }
        let symbols = encode_itinerary(&traj, &p);
        assert_eq!(symbols.len(), 100);
        assert!(symbols.iter().all(|&s| s <= 1));
        // For r=4, should see both symbols
        assert!(symbols.iter().any(|&s| s == 0));
        assert!(symbols.iter().any(|&s| s == 1));
    }

    #[test]
    fn test_shift_space_full() {
        let ss = ShiftSpace::full(2);
        assert!(ss.is_allowed(0, 0));
        assert!(ss.is_allowed(0, 1));
        assert!(ss.is_allowed(1, 0));
        assert!(ss.is_allowed(1, 1));
        assert_eq!(ss.allowed_transitions(0), 2);
    }

    #[test]
    fn test_shift_space_restricted() {
        let ss = ShiftSpace {
            n_symbols: 2,
            forbidden: vec![(0, 0)],
        };
        assert!(!ss.is_allowed(0, 0));
        assert!(ss.is_allowed(0, 1));
        assert_eq!(ss.allowed_transitions(0), 1);
        assert_eq!(ss.allowed_transitions(1), 2);
    }

    #[test]
    fn test_topological_entropy_full_shift() {
        // Full 2-shift: h_top = ln(2)
        // Use a longer pseudo-random sequence
        let mut state = 42_u64;
        let seq: Vec<usize> = (0..10000).map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize % 2
        }).collect();
        let h = topological_entropy(&seq, 8);
        // Should be close to ln(2) ≈ 0.693
        assert!((h - 2.0f64.ln()).abs() < 0.15, "h_top should be ~ln(2), got {}", h);
    }

    #[test]
    fn test_shannon_entropy() {
        // Alternating sequence: entropy = 0 (deterministic)
        let seq: Vec<usize> = (0..100).map(|i| i % 2).collect();
        let h = shannon_entropy_rate(&seq);
        // Alternating 0,1,0,1... is deterministic: entropy should be ~0
        assert!(h.abs() < 0.01, "Deterministic sequence entropy ~0, got {}", h);
    }

    #[test]
    fn test_shannon_entropy_random() {
        // Pseudo-random-ish sequence with 3 symbols using LCG
        let mut state = 123_u64;
        let seq: Vec<usize> = (0..10000).map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize % 3
        }).collect();
        let h = shannon_entropy_rate(&seq);
        // Should have entropy close to ln(3)
        assert!(h > 0.8, "Non-trivial sequence should have entropy > 0.8, got {}", h);
    }
}
