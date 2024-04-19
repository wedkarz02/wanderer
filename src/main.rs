use matrix::*;
use std::{env, error::Error, fs, process};

pub mod matrix;
pub mod monte_carlo;

fn dump_results(n: usize, eps: f64, max_iter: usize) -> Result<(), Box<dyn Error>> {
    let mat = Matrix::init_default_path(n);
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let x0 = vec![0f64; n];

    let jacobi_result = mat.jacobi(&b, &x0, eps, max_iter);
    let seidel_result = mat.gauss_seidel(&b, &x0, eps, max_iter);
    let gauss_result = mat.gaussian(&b)?;
    let gauss_partial_pivot_result = mat.gaussian_partial_pivot(&b)?;

    Matrix::vec_to_file(&jacobi_result, "dump/jacobi_result.csv")?;
    Matrix::vec_to_file(&seidel_result, "dump/seidel_result.csv")?;
    Matrix::vec_to_file(&gauss_result, "dump/gauss_result.csv")?;
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
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let x0 = vec![0f64; n];

    let jacobi_result = mat.jacobi(&b, &x0, eps, max_iter)[starting_pos];
    let seidel_result = mat.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
    let gauss_result = mat.gaussian(&b)?[starting_pos];
    let gauss_partial_pivot_result = mat.gaussian_partial_pivot(&b)?[starting_pos];

    let mc_result = monte_carlo::simulate_walk(n, starting_pos, mc_max_iter);

    let out = format!("Monte carlo: {}\nJacobi: {}\nGauss-Seidel: {}\nGauss (no pivot): {}\nGauss (partial pivot): {}\n",
        mc_result,
        jacobi_result,
        seidel_result,
        gauss_result,
        gauss_partial_pivot_result);

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
        _ => eprintln!("Unrecognised optional argument"),
    }
}
