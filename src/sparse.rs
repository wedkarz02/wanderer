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

    pub fn init_default_path(size: usize) -> Self {
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
}
