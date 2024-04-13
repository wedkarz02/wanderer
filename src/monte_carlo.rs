use rand::prelude::*;

/// Returns true if the wanderer got home safely
/// and false if he fell into the sewage well.
///
/// Panics if start is not less than the size.
pub fn walk(size: usize, start: usize) -> bool {
    assert!(start < size);

    let mut current_pos = start;
    loop {
        if current_pos == 0 {
            return true;
        } else if current_pos >= size - 1 {
            return false;
        }

        if thread_rng().gen::<bool>() {
            current_pos += 1;
        } else {
            current_pos -= 1;
        }
    }
}

pub fn simulate_walk(size: usize, start: usize, n: u32) -> f64 {
    let mut prob = 0f64;

    for _ in 0..n {
        if walk(size, start) {
            prob += 1f64;
        }
    }

    prob / n as f64
}
