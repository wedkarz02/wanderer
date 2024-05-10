use rand::prelude::*;

use crate::Config;

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

pub fn walk_park(config: &Config) -> bool {
    let mut pos = config.starting_pos + 1;

    loop {
        let mut possible_alleys = vec![];
        for alley in &config.alleys {
            if alley.a.id == pos {
                possible_alleys.push(alley.b.id);
            } else if alley.b.id == pos {
                possible_alleys.push(alley.a.id);
            }
        }

        let idx = thread_rng().gen_range(0..possible_alleys.len());
        for alley in &config.alleys {
            if alley.a.id == possible_alleys[idx] && alley.b.id == pos {
                if walk(alley.length + 1, alley.length - 1) {
                    pos = alley.a.id;
                    if alley.a.exit {
                        return true;
                    } else if alley.a.well {
                        return false;
                    }
                }
            } else if alley.b.id == possible_alleys[idx] && alley.a.id == pos {
                if walk(alley.length + 1, alley.length - 1) {
                    pos = alley.b.id;
                    if alley.b.exit {
                        return true;
                    } else if alley.b.well {
                        return false;
                    }
                }
            }
        }
    }
}

pub fn simulate_park_walk(config: &Config, max_iter: usize) -> f64 {
    let mut prob = 0f64;

    for _ in 0..max_iter {
        if walk_park(&config) {
            prob += 1f64;
        }
    }

    prob / max_iter as f64
}
