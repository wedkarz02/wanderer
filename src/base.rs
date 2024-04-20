use std::error::Error;

#[derive(Debug)]
pub enum MatrixError {
    SizeError,
    ZeroPivotError,
}

impl Error for MatrixError {}

impl std::fmt::Display for MatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SizeError => writeln!(f, "invalid matrix size"),
            Self::ZeroPivotError => writeln!(f, "zero pivot - partial pivoting is required"),
        }
    }
}

pub trait MatrixBase {
    fn init_default_path(size: usize) -> Self;
    fn to_file(&self, file_path: &'static str) -> Result<(), Box<dyn Error>>;
    fn jacobi(&self, b: &Vec<f64>, x0: &Vec<f64>, eps: f64, max_iter: usize) -> Vec<f64>;
    fn gaussian(&self, b: &Vec<f64>) -> Result<Vec<f64>, MatrixError>;
    fn partial_pivot(&self, b: &Vec<f64>) -> (Self, Vec<f64>)
    where
        Self: Sized;
    fn gaussian_partial_pivot(&self, b: &Vec<f64>) -> Result<Vec<f64>, MatrixError>;
    fn gauss_seidel(&self, b: &Vec<f64>, x0: &Vec<f64>, eps: f64, max_iter: usize) -> Vec<f64>;
}
