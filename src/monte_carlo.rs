use rand::prelude::*;

/// Returns true if the wanderer got home safely
/// and false if he fell into the sewage well.
///
/// Panics if start is not less than the size.
pub fn walk(n: usize, start: usize) -> bool {
    assert!(start < n);

    let mut current_pos = start;
    loop {
        if current_pos == 0 {
            return true;
        } else if current_pos >= n - 1 {
            return false;
        }

        if thread_rng().gen::<bool>() {
            current_pos += 1;
        } else {
            current_pos -= 1;
        }
    }
}

pub fn simulate_walk(n: usize, start: usize, max_iter: usize) -> f64 {
    let mut prob = 0f64;

    for _ in 0..max_iter {
        if walk(n, start) {
            prob += 1f64;
        }
    }

    prob / max_iter as f64
}
