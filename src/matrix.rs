pub struct Matrix {
    pub rows: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new() -> Self {
        Self { rows: vec![vec![]] }
    }

    pub fn from_size(rows: usize, cols: usize) -> Self {
        Self {
            rows: vec![vec![0.0; cols]; rows],
        }
    }
}

impl std::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in &self.rows {
            writeln!(f, "{:?}", row)?
        }
        Ok(())
    }
}
