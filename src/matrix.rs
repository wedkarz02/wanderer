use std::error::Error;
use std::{fs::File, io::Write};

use crate::base::*;

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    pub rows: Vec<Vec<f64>>,
}

impl std::fmt::Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in &self.rows {
            writeln!(f, "{:>8.4?}", row)?
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

    pub fn vec_to_file(v: &Vec<f64>, file_path: &'static str) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(file_path)?;
        let out: Vec<String> = v.iter().map(|x| x.to_string()).collect();
        let out = out.join("\n");
        file.write_all(out.as_bytes())?;
        Ok(())
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

    // https://en.wikipedia.org/wiki/Dot_product
    pub fn dot_product(x: &[f64], y: &[f64]) -> f64 {
        return x.iter().zip(y.iter()).map(|(&a, &b)| a * b).sum();
    }
}

impl MatrixBase for Matrix {
    fn init_default_path(size: usize) -> Self {
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

    fn to_file(&self, file_path: &'static str) -> Result<(), Box<dyn Error>> {
        let mut file = File::create(file_path)?;

        for row in &self.rows {
            let row_str: Vec<String> = row.iter().map(|x| x.to_string()).collect();
            let row_csv = row_str.join(";") + "\n";
            file.write_all(row_csv.as_bytes())?;
        }

        Ok(())
    }

    fn jacobi(&self, b: &Vec<f64>, x0: &Vec<f64>, eps: f64, max_iter: usize) -> Vec<f64> {
        let mut x = x0.clone();

        for it in 0..max_iter {
            let mut x_new = vec![0.0; b.len()];
            let mut error = 0f64;
            for i in 0..b.len() {
                x_new[i] = (b[i]
                    - Self::dot_product(&self.rows[i][..i], &x[..i])
                    - Self::dot_product(&self.rows[i][i + 1..], &x[i + 1..]))
                    / self.rows[i][i];

                error = error.max((x_new[i] - x[i]).abs());
            }

            if error < eps {
                println!("jacobi breaking at: {} iterations", it);
                break;
            }

            x = x_new;
        }

        x
    }

    fn gaussian(&self, b: &Vec<f64>) -> Result<Vec<f64>, MatrixError> {
        let mut a = self.clone();
        let mut b_new = b.clone();

        for i in 0..b_new.len() {
            if a.rows[i][i] == 0.0 {
                return Err(MatrixError::ZeroPivotError);
            }

            for j in (i + 1)..b_new.len() {
                let factor = a.rows[j][i] / a.rows[i][i];
                for k in i..b_new.len() {
                    a.rows[j][k] -= factor * a.rows[i][k];
                }
                b_new[j] -= factor * b_new[i];
            }
        }

        let mut out = vec![0f64; b_new.len()];
        for i in (0..b_new.len()).rev() {
            out[i] = b_new[i];
            for j in (i + 1)..b_new.len() {
                out[i] -= a.rows[i][j] * out[j];
            }
            out[i] /= a.rows[i][i];
        }

        Ok(out)
    }

    fn partial_pivot(&self, b: &Vec<f64>) -> (Self, Vec<f64>) {
        let mut a = self.clone();
        let mut b_new = b.clone();

        for i in 0..b_new.len() {
            let mut max_row = i;
            for k in (i + 1)..b_new.len() {
                if a.rows[k][i].abs() > a.rows[max_row][i].abs() {
                    max_row = k;
                }
            }

            if max_row != i {
                a.rows.swap(i, max_row);
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

        for it in 0..max_iter {
            let mut x_new = vec![0f64; b.len()];
            let mut error = 0f64;
            for i in 0..b.len() {
                x_new[i] = (b[i]
                    - Self::dot_product(&self.rows[i][..i], &x_new[..i])
                    - Self::dot_product(&self.rows[i][i + 1..], &x[i + 1..]))
                    / self.rows[i][i];

                error = error.max((x_new[i] - x[i]).abs());
            }

            if error < eps {
                println!("seidel breaking at: {} iterations", it);
                break;
            }

            x = x_new;
        }

        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply() {
        let a = Matrix::from_vecs(vec![
            vec![1.0, 2.0, 1.0],
            vec![2.0, 1.0, 1.0],
            vec![1.0, 3.0, 1.0],
        ]);
        let b = Matrix::from_vecs(vec![vec![5.0], vec![2.0], vec![4.0]]);
        let expected_product = Matrix::from_vecs(vec![vec![13.0], vec![16.0], vec![15.0]]);
        assert_eq!(a.multiply(&b).unwrap(), expected_product);
    }
}
