use crate::base::*;
use crate::Config;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Sparse {
    data: HashMap<(usize, usize), f64>,
}

impl std::fmt::Display for Sparse {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (pos, val) in self.data.iter() {
            writeln!(f, "{:?}: {}", pos, val)?
        }
        Ok(())
    }
}

impl Sparse {
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    pub fn from_size(size: usize) -> Self {
        Self {
            data: HashMap::with_capacity(size),
        }
    }

    pub fn from_vecs(vec_matrix: Vec<Vec<f64>>) -> Self {
        let mut sparse = Self::from_size(vec_matrix.len());

        for i in 0..vec_matrix.len() {
            for j in 0..vec_matrix[i].len() {
                if vec_matrix[i][j] != 0f64 {
                    sparse.data.insert((i, j), vec_matrix[i][j]);
                }
            }
        }

        sparse
    }

    pub fn from_config(cfg: &Config) -> (Self, Vec<f64>) {
        let n = cfg.inters.len();
        let mut out = Self::from_size(n);
        let mut b = vec![0f64; n];

        for i in 0..n {
            out.data.insert((i, i), 1f64);
            if cfg.inters[i].exit {
                b[i] = 1f64;
            }
        }

        for i in 0..n {
            if cfg.inters[i].exit || cfg.inters[i].well {
                continue;
            }

            let mut alley_ctr = 0f64;
            for alley in &cfg.alleys {
                if cfg.inters[i].id == alley.a.id || cfg.inters[i].id == alley.b.id {
                    alley_ctr += 1f64;
                }
            }

            for alley in &cfg.alleys {
                if cfg.inters[i].id == alley.a.id {
                    let value = if !cfg.deadends.contains(&alley.a.id) {
                        let mut tmp = -alley.get_propability() / alley_ctr;
                        if alley.b.trashcan {
                            tmp /= 2f64;
                        }
                        tmp
                    } else {
                        1f64
                    };
                    out.data.insert((i, alley.b.id - 1), value);
                } else if cfg.inters[i].id == alley.b.id {
                    let value = if !cfg.deadends.contains(&alley.b.id) {
                        let mut tmp = -alley.get_propability() / alley_ctr;
                        if alley.a.trashcan {
                            tmp /= 2f64;
                        }
                        tmp
                    } else {
                        1f64
                    };
                    out.data.insert((i, alley.a.id - 1), value);
                }
            }
        }

        (out, b)
    }

    pub fn get_value(&self, i: usize, j: usize) -> f64 {
        *self.data.get(&(i, j)).unwrap_or(&0f64)
    }
}

impl MatrixBase for Sparse {
    fn init_default_path(size: usize) -> Self {
        let mut sparse = Self::from_size(size);

        sparse.data.insert((0, 0), 1f64);
        sparse.data.insert((size - 1, size - 1), 1f64);

        for i in 1..size - 1 {
            sparse.data.insert((i, i - 1), -0.5);
            sparse.data.insert((i, i), 1f64);
            sparse.data.insert((i, i + 1), -0.5);
        }

        sparse
    }

    // This function is heavily inspired by demo C++ sparse implementation
    // written by my lecturer dr. Łukasz Kuszner.
    fn jacobi(&self, b: &Vec<f64>, x0: &Vec<f64>, eps: f64, max_iter: usize) -> Vec<f64> {
        let mut x = x0.clone();

        for it in 0..max_iter {
            let mut error = 0f64;
            let x_temp = x.clone();
            for i in 0..b.len() {
                x[i] = b[i];
            }

            for (pos, val) in &self.data {
                if pos.0 != pos.1 {
                    x[pos.0] -= x_temp[pos.1] * val;
                }
            }

            for (pos, val) in &self.data {
                if pos.0 == pos.1 {
                    x[pos.0] /= val;
                }

                error = error.max((x_temp[pos.0] - x[pos.0]).abs());
            }

            if error < eps {
                println!("sparse jacobi breaking at: {} iterations", it);
                break;
            }
        }

        x
    }

    fn gaussian(&self, b: &Vec<f64>) -> Result<Vec<f64>, MatrixError> {
        let mut a = self.clone();
        let mut b_new = b.clone();

        for i in 0..b_new.len() {
            let pivot = match a.data.get(&(i, i)) {
                Some(&val) => val,
                None => return Err(MatrixError::ZeroPivotError),
            };

            for j in (i + 1)..b_new.len() {
                let factor = match a.data.get(&(j, i)) {
                    Some(&val) => val / pivot,
                    None => continue,
                };

                for k in i..b_new.len() {
                    let entry = a.get_value(i, k);
                    if let Some(val) = a.data.get_mut(&(j, k)) {
                        *val -= factor * entry;
                    }
                }
                b_new[j] -= factor * b_new[i];
            }
        }

        let mut out = vec![0f64; b_new.len()];
        for i in (0..b_new.len()).rev() {
            out[i] = b_new[i];
            for j in (i + 1)..b_new.len() {
                out[i] -= a.get_value(i, j) * out[j];
            }
            out[i] /= a.get_value(i, i);
            if out[i].is_nan() {
                return Err(MatrixError::Unsolvable);
            }
        }

        Ok(out)
    }

    fn partial_pivot(&self, b: &Vec<f64>) -> (Self, Vec<f64>) {
        let mut a = self.clone();
        let mut b_new = b.clone();

        for i in 0..b_new.len() {
            let mut max_row = i;
            for k in (i + 1)..b_new.len() {
                match (a.data.get(&(k, i)), a.data.get(&(max_row, i))) {
                    (Some(&val), Some(&max_val)) => {
                        if val.abs() > max_val.abs() {
                            max_row = k;
                        }
                    }
                    (Some(_), None) => max_row = k,
                    _ => {}
                }
            }

            if max_row != i {
                for j in 0..b_new.len() {
                    if let Some(&val) = a.data.get(&(i, j)) {
                        a.data.insert((i, j), a.get_value(max_row, j));
                        a.data.insert((max_row, j), val);
                    } else if let Some(&max_val) = a.data.get(&(max_row, j)) {
                        a.data.insert((i, j), max_val);
                    }
                }
                b_new.swap(i, max_row);
            }
        }

        (a, b_new)
    }

    fn gaussian_partial_pivot(&self, b: &Vec<f64>) -> Result<Vec<f64>, MatrixError> {
        let (a, b_new) = self.partial_pivot(b);
        return a.gaussian(&b_new);
    }

    fn gauss_seidel(&self, b: &Vec<f64>, x0: &Vec<f64>, eps: f64, max_iter: usize) -> Vec<f64> {
        let mut x = x0.clone();

        // FIXME: This is ridiculously slow for some reason
        //        590 seconds for 1k x 1k matrix, 1e-6 eps and 1mil max_iter

        for it in 0..max_iter {
            let mut x_new = vec![0.0; b.len()];
            let mut error = 0f64;

            for i in 0..b.len() {
                for (pos, val) in &self.data {
                    if pos.0 != pos.1 && pos.0 == i {
                        x_new[pos.0] += val * if pos.1 < i { x_new[pos.1] } else { x[pos.1] };
                    }
                }

                x_new[i] = (b[i] - x_new[i]) / self.get_value(i, i);
                error = error.max((x_new[i] - x[i]).abs());
            }

            if error < eps {
                println!("Gauss-Seidel breaking at: {} iterations", it);
                break;
            }

            x = x_new;
        }

        x
    }
}
