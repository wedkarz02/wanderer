use base::*;
use matrix::*;
use sparse::*;
use std::time::Instant;
use std::{env, error::Error, fs, process};

pub mod base;
pub mod matrix;
pub mod monte_carlo;
pub mod sparse;

fn dump_results(n: usize, eps: f64, max_iter: usize) -> Result<(), Box<dyn Error>> {
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

    Matrix::vec_to_file(&jacobi_result, "dump/jacobi_result.csv")?;
    Matrix::vec_to_file(&jacobi_sparse_result, "dump/jacobi_sparse_result.csv")?;
    Matrix::vec_to_file(&seidel_result, "dump/seidel_result.csv")?;
    Matrix::vec_to_file(&gauss_result, "dump/gauss_result.csv")?;
    Matrix::vec_to_file(&gauss_sparse_result, "dump/gauss_sparse_result.csv")?;
    Matrix::vec_to_file(
        &gauss_partial_pivot_result,
        "dump/gauss_partial_pivot_result.csv",
    )?;

    Ok(())
}

fn mc_compare(
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

    let gauss_start = Instant::now();
    let gauss_result = mat.gaussian(&b)?[starting_pos];
    let gauss_elapsed = gauss_start.elapsed();

    let g_sparse_start = Instant::now();
    let g_sparse_result = sparse_mat.gaussian(&b)?[starting_pos];
    let g_sparse_elapsed = g_sparse_start.elapsed();

    let gpp_start = Instant::now();
    let gauss_partial_pivot_result = mat.gaussian_partial_pivot(&b)?[starting_pos];
    let gpp_elapsed = gpp_start.elapsed();

    let mc_start = Instant::now();
    let mc_result = monte_carlo::simulate_walk(n, starting_pos, mc_max_iter);
    let mc_elapsed = mc_start.elapsed();

    let out = format!(
        "Monte carlo: {} in {:?}\n\
        Jacobi: {} in {:?}\n\
        Sparse Jacobi: {} in {:?}\n\
        Gauss-Seidel: {} in {:?}\n\
        Gauss (no pivot): {} in {:?}\n\
        Sparse Gauss (no pivot): {} in {:?}\n\
        Gauss (partial pivot): {} in {:?}\n",
        mc_result,
        mc_elapsed,
        jacobi_result,
        jacobi_elapsed,
        jacobi_sparse_result,
        jacobi_sparse_elapsed,
        seidel_result,
        seidel_elapsed,
        gauss_result,
        gauss_elapsed,
        g_sparse_result,
        g_sparse_elapsed,
        gauss_partial_pivot_result,
        gpp_elapsed,
    );

    fs::write("dump/result_comparison.txt", out)?;

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Nothing to do");
        process::exit(0);
    }

    let n = 1000;
    let eps = 1e-12;
    let max_iter = 100_000;
    let mc_max_iter = 10_000;
    let starting_pos = n / 2;

    match args[1].as_str() {
        "dump" => {
            if let Err(e) = dump_results(n, eps, max_iter) {
                eprintln!("Test data dump failed: {}", e);
            }
        }
        "compare" => {
            if let Err(e) = mc_compare(n, eps, max_iter, mc_max_iter, starting_pos) {
                eprintln!("Comparison dump failed: {}", e);
            }
        }
        "sparse" => {
            let n = 1000;
            let mut b = vec![0f64; n];
            b[0] = 1f64;

            let mat: Matrix = MatrixBase::init_default_path(n);
            let start = Instant::now();
            let res = mat.gaussian(&b);
            let time = start.elapsed();

            println!("{:?}\n{:#?}", time, res);
        }
        _ => eprintln!("Unrecognised optional argument"),
    }
}
