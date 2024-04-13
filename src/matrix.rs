#[derive(Debug)]
pub enum MatrixError {
    SizeError,
}

impl std::fmt::Display for MatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SizeError => writeln!(f, "invalid matrix size"),
        }
    }
}

pub struct Matrix {
    pub rows: Vec<Vec<f64>>,
}

impl std::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in &self.rows {
            writeln!(f, "{:4?}", row)?
        }
        Ok(())
    }
}

impl Matrix {
    pub fn new() -> Self {
        Self { rows: vec![] }
    }

    pub fn from_size(rows: usize, cols: usize) -> Self {
        Self {
            rows: vec![vec![0f64; cols]; rows],
        }
    }

    pub fn from_vecs(vec_matrix: Vec<Vec<f64>>) -> Self {
        Self { rows: vec_matrix }
    }

    pub fn init_default_path(size: usize) -> Self {
        let mut out = Self::from_size(size, size);

        out.rows[0][0] = 1f64;
        out.rows[size - 1][size - 1] = 1f64;

        for i in 1..size - 1 {
            out.rows[i][i - 1] = -0.5;
            out.rows[i][i] = 1f64;
            out.rows[i][i + 1] = -0.5;
        }

        out
    }

    pub fn multiply(&self, other: &Self) -> Result<Self, MatrixError> {
        if self.rows[0].len() != other.rows.len() {
            return Err(MatrixError::SizeError);
        }

        let mut out = Self::from_size(self.rows.len(), other.rows[0].len());

        for i in 0..self.rows.len() {
            for j in 0..other.rows[0].len() {
                for k in 0..self.rows[0].len() {
                    out.rows[i][j] += self.rows[i][k] * other.rows[k][j];
                }
            }
        }

        Ok(out)
    }

    pub fn dot_product(x: &[f64], y: &[f64]) -> f64 {
        return x.iter().zip(y.iter()).map(|(&a, &b)| a * b).sum();
    }

    // https://en.wikipedia.org/wiki/Jacobi_method
    pub fn jacobi(&self, b: &Vec<f64>, x0: &Vec<f64>, max_iter: usize) -> Vec<f64> {
        let mut x = x0.to_vec();
        let mut x_new = vec![0.0; self.rows.len()];

        for _ in 0..max_iter {
            for i in 0..self.rows.len() {
                x_new[i] = (b[i]
                    - Self::dot_product(&self.rows[i][..i], &x[..i])
                    - Self::dot_product(&self.rows[i][i + 1..], &x[i + 1..]))
                    / self.rows[i][i];
            }
            x = x_new.clone();
        }

        x_new
    }
}
