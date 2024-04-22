use crate::base::*;
use crate::matrix::*;
use crate::sparse::*;

use std::{error::Error, fs, time::Instant};

fn compare_default(
    n: usize,
    eps: f64,
    max_iter: usize,
    starting_pos: usize,
) -> Result<(String, String, String, String), Box<dyn Error>> {
    let mat = Matrix::init_default_path(n);
    let sparse_mat = Sparse::init_default_path(n);
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let x0 = vec![0f64; n];

    let jacobi_start = Instant::now();
    let jacobi_result = mat.jacobi(&b, &x0, eps, max_iter)[starting_pos];
    let jacobi_elapsed = jacobi_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(jacobi_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let jacobi_sparse_start = Instant::now();
    let jacobi_sparse_result = sparse_mat.jacobi(&b, &x0, eps, max_iter)[starting_pos];
    let jacobi_sparse_elapsed = jacobi_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(jacobi_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let seidel_start = Instant::now();
    let seidel_result = mat.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
    let seidel_elapsed = seidel_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(seidel_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_start = Instant::now();
    let gauss_result = mat.gaussian(&b)?[starting_pos];
    let gauss_elapsed = gauss_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gauss_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let g_sparse_start = Instant::now();
    let g_sparse_result = sparse_mat.gaussian(&b)?[starting_pos];
    let g_sparse_elapsed = g_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(g_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gpp_start = Instant::now();
    let gauss_partial_pivot_result = mat.gaussian_partial_pivot(&b)?[starting_pos];
    let gpp_elapsed = gpp_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gpp_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gpp_sparse_start = Instant::now();
    let gpp_sparse_result = sparse_mat.gaussian_partial_pivot(&b)?[starting_pos];
    let gpp_sparse_elapsed = gpp_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gpp_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let non_sparse_row = format!(
        "{:.6?};{:.6?};{:.6?};{:.6?};{:.6?}",
        n, jacobi_elapsed, seidel_elapsed, gauss_elapsed, gpp_elapsed,
    );

    let non_sparse_results = format!(
        "{};{};{};{};{}",
        n, jacobi_result, seidel_result, gauss_result, gauss_partial_pivot_result,
    );

    let sparse_row = format!(
        "{:.6?};{:.6?};{:.6?};{:.6?}",
        n, jacobi_sparse_elapsed, g_sparse_elapsed, gpp_sparse_elapsed,
    );

    let sparse_results = format!(
        "{};{};{};{}",
        n, jacobi_sparse_result, g_sparse_result, gpp_sparse_result
    );

    Ok((
        non_sparse_row,
        non_sparse_results,
        sparse_row,
        sparse_results,
    ))
}

pub fn incremental_compare() {
    let mut ns_row_lines = Vec::new();
    let mut ns_res_lines = Vec::new();
    let mut s_row_lines = Vec::new();
    let mut s_res_lines = Vec::new();

    ns_row_lines.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    ns_res_lines.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    s_row_lines.push(String::from("n;jacobi;gauss;gauss_pivot"));
    s_res_lines.push(String::from("n;jacobi;gauss;gauss_pivot"));

    for n in (10..=300).step_by(10) {
        let starting_pos = n / 2;
        let eps = 1e-16;
        let max_iter = 1_000_000;

        let (ns_row, ns_res, s_row, s_res) =
            compare_default(n, eps, max_iter, starting_pos).unwrap();

        ns_row_lines.push(ns_row);
        ns_res_lines.push(ns_res);
        s_row_lines.push(s_row);
        s_res_lines.push(s_res);
    }

    let ns_row_str = ns_row_lines.join("\n");
    let ns_res_str = ns_res_lines.join("\n");
    let s_row_str = s_row_lines.join("\n");
    let s_res_str = s_res_lines.join("\n");

    fs::write("dump/default_no_sparse_time.csv", ns_row_str).unwrap();
    fs::write("dump/default_no_sparse_results.csv", ns_res_str).unwrap();
    fs::write("dump/default_sparse_time.csv", s_row_str).unwrap();
    fs::write("dump/default_sparse_results.csv", s_res_str).unwrap();
}
