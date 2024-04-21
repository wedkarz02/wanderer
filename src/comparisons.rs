use crate::base::*;
use crate::matrix::*;
use crate::monte_carlo;
use crate::sparse::*;

use std::{error::Error, fs, time::Instant};

pub fn dump_results(n: usize, eps: f64, max_iter: usize) -> Result<(), Box<dyn Error>> {
    let mat: Matrix = MatrixBase::init_default_path(n);
    let sparse_mat: Sparse = Sparse::init_default_path(n);
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let x0 = vec![0f64; n];

    let jacobi_result = mat.jacobi(&b, &x0, eps, max_iter);
    let jacobi_sparse_result = sparse_mat.jacobi(&b, &x0, eps, max_iter);
    let seidel_result = mat.gauss_seidel(&b, &x0, eps, max_iter);
    let gauss_result = mat.gaussian(&b)?;
    let gauss_sparse_result = sparse_mat.gaussian(&b)?;
    let gauss_partial_pivot_result = mat.gaussian_partial_pivot(&b)?;
    let gpp_sparse_result = sparse_mat.gaussian_partial_pivot(&b)?;

    Matrix::vec_to_file(&jacobi_result, "dump/jacobi_result.csv")?;
    Matrix::vec_to_file(&jacobi_sparse_result, "dump/jacobi_sparse_result.csv")?;
    Matrix::vec_to_file(&seidel_result, "dump/seidel_result.csv")?;
    Matrix::vec_to_file(&gauss_result, "dump/gauss_result.csv")?;
    Matrix::vec_to_file(&gauss_sparse_result, "dump/gauss_sparse_result.csv")?;
    Matrix::vec_to_file(
        &gauss_partial_pivot_result,
        "dump/gauss_partial_pivot_result.csv",
    )?;
    Matrix::vec_to_file(&gpp_sparse_result, "dump/gauss_pp_sparse_result.csv")?;

    Ok(())
}

pub fn mc_compare(
    n: usize,
    eps: f64,
    max_iter: usize,
    mc_max_iter: usize,
    starting_pos: usize,
) -> Result<(), Box<dyn Error>> {
    let mat = Matrix::init_default_path(n);
    let sparse_mat = Sparse::init_default_path(n);
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let x0 = vec![0f64; n];

    let jacobi_start = Instant::now();
    let jacobi_result = mat.jacobi(&b, &x0, eps, max_iter)[starting_pos];
    let jacobi_elapsed = jacobi_start.elapsed();

    let jacobi_sparse_start = Instant::now();
    let jacobi_sparse_result = sparse_mat.jacobi(&b, &x0, eps, max_iter)[starting_pos];
    let jacobi_sparse_elapsed = jacobi_sparse_start.elapsed();

    let seidel_start = Instant::now();
    let seidel_result = mat.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
    let seidel_elapsed = seidel_start.elapsed();

    let s_sparse_start = Instant::now();
    let s_sparse_result = sparse_mat.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
    let s_sparse_elapsed = s_sparse_start.elapsed();

    let gauss_start = Instant::now();
    let gauss_result = mat.gaussian(&b)?[starting_pos];
    let gauss_elapsed = gauss_start.elapsed();

    let g_sparse_start = Instant::now();
    let g_sparse_result = sparse_mat.gaussian(&b)?[starting_pos];
    let g_sparse_elapsed = g_sparse_start.elapsed();

    let gpp_start = Instant::now();
    let gauss_partial_pivot_result = mat.gaussian_partial_pivot(&b)?[starting_pos];
    let gpp_elapsed = gpp_start.elapsed();

    let gpp_sparse_start = Instant::now();
    let gpp_sparse_result = sparse_mat.gaussian_partial_pivot(&b)?[starting_pos];
    let gpp_sparse_elapsed = gpp_sparse_start.elapsed();

    let mc_start = Instant::now();
    let mc_result = monte_carlo::simulate_walk(n, starting_pos, mc_max_iter);
    let mc_elapsed = mc_start.elapsed();

    let out = format!(
        "Monte carlo: {} in {:?}\n\
        Jacobi: {} in {:?}\n\
        Sparse Jacobi: {} in {:?}\n\
        Gauss-Seidel: {} in {:?}\n\
        Sparse Gauss-Seidel: {} in {:?}\n\
        Gauss (no pivot): {} in {:?}\n\
        Sparse Gauss (no pivot): {} in {:?}\n\
        Gauss (partial pivot): {} in {:?}\n\
        Sparse Gauss (partial pivot): {} in {:?}",
        mc_result,
        mc_elapsed,
        jacobi_result,
        jacobi_elapsed,
        jacobi_sparse_result,
        jacobi_sparse_elapsed,
        seidel_result,
        seidel_elapsed,
        s_sparse_result,
        s_sparse_elapsed,
        gauss_result,
        gauss_elapsed,
        g_sparse_result,
        g_sparse_elapsed,
        gauss_partial_pivot_result,
        gpp_elapsed,
        gpp_sparse_result,
        gpp_sparse_elapsed
    );

    fs::write("dump/result_comparison.txt", out)?;

    Ok(())
}
