#[derive(Debug)]
pub struct SizeError;

impl std::fmt::Display for SizeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "invalid matrix size")
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
        Self { rows: vec![vec![]] }
    }

    pub fn from_size(rows: usize, cols: usize) -> Self {
        Self {
            rows: vec![vec![0f64; cols]; rows],
        }
    }

    pub fn multiply(&self, other: &Self) -> Result<Self, SizeError> {
        if self.rows[0].len() != other.rows.len() {
            return Err(SizeError);
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
}
