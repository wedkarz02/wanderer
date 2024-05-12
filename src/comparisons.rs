use crate::base::*;
use crate::matrix::*;
use crate::monte_carlo;
use crate::sparse::*;
use crate::Config;

use std::process;
use std::{error::Error, fs, time::Instant};

pub fn compare_vecs(a: &Vec<f64>, b: &Vec<f64>, eps: f64) -> bool {
    if a.len() != b.len() {
        return false;
    }

    for i in 0..a.len() {
        if (a[i] - b[i]).abs() > eps {
            return false;
        }
    }

    true
}

pub fn compare_config(
    cfg: &Config,
    eps: f64,
    max_iter: usize,
) -> Result<(String, String, String, String), MatrixError> {
    let (mat, b) = Matrix::from_config(cfg);
    let (sparse, _) = Sparse::from_config(cfg);
    let x0 = vec![0f64; b.len()];

    let gpp_start = Instant::now();
    let gpp_result = mat.gaussian_partial_pivot(&b);
    let gpp_elapsed = gpp_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gpp_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gpp_result = match gpp_result {
        Ok(values) => values,
        Err(e) => return Err(e),
    };

    let gpp_sparse_start = Instant::now();
    let gpp_sparse_res = sparse.gaussian_partial_pivot(&b);
    let gpp_sparse_elapsed = gpp_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gpp_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gpp_sparse_res = match gpp_sparse_res {
        Ok(values) => values,
        Err(e) => return Err(e),
    };

    let jacobi_start = Instant::now();
    let jacobi_res = mat.jacobi(&b, &x0, eps, max_iter);
    let jacobi_elapsed = jacobi_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(jacobi_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let jacobi_sparse_start = Instant::now();
    let jacobi_sparse_res = sparse.jacobi(&b, &x0, eps, max_iter);
    let jacobi_sparse_elapsed = jacobi_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(jacobi_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let seidel_start = Instant::now();
    let seidel_res = mat.gauss_seidel(&b, &x0, eps, max_iter);
    let seidel_elapsed = seidel_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(seidel_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_start = Instant::now();
    let gauss_res = mat.gaussian(&b);
    let gauss_elapsed = gauss_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gauss_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_res = match gauss_res {
        Ok(values) => values,
        Err(e) => return Err(e),
    };

    let gauss_sparse_start = Instant::now();
    let gauss_sparse_res = sparse.gaussian(&b);
    let gauss_sparse_elapsed = gauss_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gauss_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_sparse_res = match gauss_sparse_res {
        Ok(values) => values,
        Err(e) => return Err(e),
    };

    let seidel_sparse_start = Instant::now();
    let seidel_sparse_res = mat.gauss_seidel(&b, &x0, eps, max_iter);
    let seidel_sparse_elapsed = seidel_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(seidel_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let res_ns_line = format!(
        "{};{};{};{};{}",
        b.len(),
        jacobi_res[cfg.starting_pos],
        seidel_res[cfg.starting_pos],
        gauss_res[cfg.starting_pos],
        gpp_result[cfg.starting_pos],
    );
    let time_ns_line = format!(
        "{};{};{};{};{}",
        b.len(),
        jacobi_elapsed,
        seidel_elapsed,
        gauss_elapsed,
        gpp_elapsed
    );

    let res_s_line = format!(
        "{};{};{};{};{}",
        b.len(),
        jacobi_sparse_res[cfg.starting_pos],
        seidel_sparse_res[cfg.starting_pos],
        gauss_sparse_res[cfg.starting_pos],
        gpp_sparse_res[cfg.starting_pos]
    );
    let time_s_line = format!(
        "{};{};{};{};{}",
        b.len(),
        jacobi_sparse_elapsed,
        seidel_sparse_elapsed,
        gauss_sparse_elapsed,
        gpp_sparse_elapsed
    );

    Ok((res_ns_line, time_ns_line, res_s_line, time_s_line))
}

pub fn incremental_compare_config(
    n: usize,
    step: usize,
    cfg: Option<&Config>,
) -> Result<(), Box<dyn Error>> {
    let mut results_ns = Vec::new();
    let mut times_ns = Vec::new();
    let mut results_s = Vec::new();
    let mut times_s = Vec::new();

    results_ns.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    times_ns.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    results_s.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    times_s.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));

    let eps = 1e-16;
    let max_iter = 10_000;

    let mut res_ns_line;
    let mut time_ns_line;
    let mut res_s_line;
    let mut time_s_line;

    for i in (step..=n * step).step_by(step) {
        loop {
            let config = match cfg {
                Some(c) => c.clone(),
                None => {
                    let mut tmp_gauss: Result<Vec<f64>, MatrixError> = Err(MatrixError::Unsolvable);
                    while let Err(_) = tmp_gauss {
                        crate::gen_config(i, (3 * i) / 2)?;
                        let sets = crate::parse_config("tmp.config");
                        let config = Config::build(sets);
                        let (sparse, b) = Sparse::from_config(&config);
                        tmp_gauss = sparse.gaussian(&b);
                    }
                    let sets = crate::parse_config("tmp.config");
                    Config::build(sets).clone()
                }
            };

            (res_ns_line, time_ns_line, res_s_line, time_s_line) =
                compare_config(&config, eps, max_iter)?;

            if !res_ns_line.contains("-") || !res_s_line.contains("-") {
                break;
            } else {
                println!("{:#?}", config);
            }
        }

        results_ns.push(res_ns_line);
        times_ns.push(time_ns_line);
        results_s.push(res_s_line);
        times_s.push(time_s_line);
    }

    let results_ns_str = results_ns.join("\n");
    let times_ns_str = times_ns.join("\n");
    let results_s_str = results_s.join("\n");
    let times_s_str = times_s.join("\n");

    match cfg {
        Some(_) => {
            fs::write("dump/single_config_res_no_sparse.csv", results_ns_str)?;
            fs::write("dump/single_config_time_no_sparse.csv", times_ns_str)?;
            fs::write("dump/single_config_res_sparse.csv", results_s_str)?;
            fs::write("dump/single_config_time_sparse.csv", times_s_str)?;
        }
        None => {
            fs::write("dump/config_res_no_sparse.csv", results_ns_str)?;
            fs::write("dump/config_time_no_sparse.csv", times_ns_str)?;
            fs::write("dump/config_res_sparse.csv", results_s_str)?;
            fs::write("dump/config_time_sparse.csv", times_s_str)?;
        }
    }

    Ok(())
}

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

    let seidel_sparse_start = Instant::now();
    let seidel_sparse_result = mat.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
    let seidel_sparse_elapsed = seidel_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(seidel_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let non_sparse_row = format!(
        "{};{};{};{};{}",
        n, jacobi_elapsed, seidel_elapsed, gauss_elapsed, gpp_elapsed,
    );

    let non_sparse_results = format!(
        "{};{};{};{};{}",
        n, jacobi_result, seidel_result, gauss_result, gauss_partial_pivot_result,
    );

    let sparse_row = format!(
        "{};{};{};{};{}",
        n, jacobi_sparse_elapsed, seidel_sparse_elapsed, g_sparse_elapsed, gpp_sparse_elapsed,
    );

    let sparse_results = format!(
        "{};{};{};{};{}",
        n, jacobi_sparse_result, seidel_sparse_result, g_sparse_result, gpp_sparse_result
    );

    Ok((
        non_sparse_row,
        non_sparse_results,
        sparse_row,
        sparse_results,
    ))
}

pub fn incremental_compare_default() {
    let mut ns_row_lines = Vec::new();
    let mut ns_res_lines = Vec::new();
    let mut s_row_lines = Vec::new();
    let mut s_res_lines = Vec::new();

    ns_row_lines.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    ns_res_lines.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    s_row_lines.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));
    s_res_lines.push(String::from("n;jacobi;seidel;gauss;gauss_pivot"));

    let eps = 1e-16;
    let max_iter = 1_000;

    for n in (10..=300).step_by(10) {
        let starting_pos = n / 2;

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

pub fn incremental_verify_mc(n: usize, cfg: Option<&Config>) -> Result<(), Box<dyn Error>> {
    let mut lines = Vec::new();
    lines.push(String::from("n;gauss_pivot;mc;err"));

    for i in (10..=n * 10).step_by(10) {
        let config = match cfg {
            Some(c) => c.clone(),
            None => {
                let mut tmp_res = 0f64;
                while tmp_res == 0f64 || tmp_res == 1f64 {
                    crate::gen_config(i, (3 * i) / 2)?;
                    let sets = crate::parse_config("tmp.config");
                    let config = Config::build(sets);
                    let (sparse, b) = Sparse::from_config(&config);
                    tmp_res = sparse.gaussian_partial_pivot(&b).unwrap()[config.starting_pos];
                }
                let sets = crate::parse_config("tmp.config");
                Config::build(sets).clone()
            }
        };

        let (sparse, b) = Sparse::from_config(&config);
        let sparse_res = match sparse.gaussian_partial_pivot(&b) {
            Ok(values) => values[config.starting_pos],
            Err(e) => return Err(Box::new(e)),
        };

        let mc_res = monte_carlo::simulate_park_walk(&config, i * 200);

        lines.push(format!(
            "{};{};{};{}",
            i,
            sparse_res,
            mc_res,
            (sparse_res - mc_res).abs(),
        ));
    }

    let lines_str = lines.join("\n");
    fs::write("dump/mc_verify.csv", lines_str).unwrap();

    Ok(())
}

pub fn check_results(config: &Config) {
    let (mat, b) = Matrix::from_config(&config);
    let (sparse, _) = Sparse::from_config(&config);
    let x0 = vec![0f64; b.len()];
    let eps = 1e-16;
    let max_iter = 10_000;

    let gpp_res = match mat.gaussian_partial_pivot(&b) {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };
    let gpp_sub = mat.multiply_by_vec(&gpp_res).unwrap();
    let gpp_sparse_res = match sparse.gaussian_partial_pivot(&b) {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };
    let gpp_sparse_sub = mat.multiply_by_vec(&gpp_sparse_res).unwrap();
    let jacobi_res = mat.jacobi(&b, &x0, eps, max_iter);
    let jacobi_sub = mat.multiply_by_vec(&jacobi_res).unwrap();
    let jacobi_sparse_res = sparse.jacobi(&b, &x0, eps, max_iter);
    let jacobi_sparse_sub = mat.multiply_by_vec(&jacobi_sparse_res).unwrap();
    let seidel_res = mat.gauss_seidel(&b, &x0, eps, max_iter);
    let seidel_sub = mat.multiply_by_vec(&seidel_res).unwrap();
    let seidel_sparse_res = sparse.gauss_seidel(&b, &x0, eps, max_iter);
    let seidel_sparse_sub = mat.multiply_by_vec(&seidel_sparse_res).unwrap();
    let gauss_res = match mat.gaussian(&b) {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };
    let gauss_sub = mat.multiply_by_vec(&gauss_res).unwrap();
    let gauss_sparse_res = match sparse.gaussian(&b) {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };
    let gauss_sparse_sub = mat.multiply_by_vec(&gauss_sparse_res).unwrap();

    if compare_vecs(&gpp_sub, &b, eps) {
        println!("\nGauss (partial pivot): Success")
    } else {
        println!("\nGauss (partial pivot): Failure")
    }
    if compare_vecs(&gpp_sparse_sub, &b, eps) {
        println!("Sparse Gauss (partial pivot): Success")
    } else {
        println!("Sparse Gauss (partial pivot): Failure")
    }
    if compare_vecs(&jacobi_sub, &b, eps) {
        println!("Jacobi: Success")
    } else {
        println!("Jacobi: Failure")
    }
    if compare_vecs(&jacobi_sparse_sub, &b, eps) {
        println!("Sparse Jacobi: Success")
    } else {
        println!("Sparse Jacobi: Failure")
    }
    if compare_vecs(&seidel_sub, &b, eps) {
        println!("Gauss-Seidel: Success")
    } else {
        println!("Gauss-Seidel: Failure")
    }
    if compare_vecs(&seidel_sparse_sub, &b, eps) {
        println!("Sparse Gauss-Seidel: Success")
    } else {
        println!("Sparse Gauss-Seidel: Failure")
    }
    if compare_vecs(&gauss_sub, &b, eps) {
        println!("Gauss: Success")
    } else {
        println!("Gauss: Failure")
    }
    if compare_vecs(&gauss_sparse_sub, &b, eps) {
        println!("Sparse Gauss: Success")
    } else {
        println!("Sparse Gauss: Failure")
    }
}

pub fn time_all(config: &Config) {
    let (mat, _) = Matrix::from_config(&config);
    let (sparse, b) = Sparse::from_config(&config);
    let x0 = vec![0f64; b.len()];
    let eps = 1e-16;
    let max_iter = 1_000;

    let gpp_start = Instant::now();
    let gpp_result = mat.gaussian_partial_pivot(&b);
    let gpp_elapsed = gpp_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gpp_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gpp_result = match gpp_result {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };

    let gpp_sparse_start = Instant::now();
    let gpp_sparse_res = sparse.gaussian_partial_pivot(&b);
    let gpp_sparse_elapsed = gpp_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gpp_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gpp_sparse_res = match gpp_sparse_res {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };

    let jacobi_start = Instant::now();
    let jacobi_res = mat.jacobi(&b, &x0, eps, max_iter);
    let jacobi_elapsed = jacobi_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(jacobi_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let jacobi_sparse_start = Instant::now();
    let jacobi_sparse_res = sparse.jacobi(&b, &x0, eps, max_iter);
    let jacobi_sparse_elapsed = jacobi_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(jacobi_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let seidel_start = Instant::now();
    let seidel_res = mat.gauss_seidel(&b, &x0, eps, max_iter);
    let seidel_elapsed = seidel_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(seidel_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let seidel_sparse_start = Instant::now();
    let seidel_sparse_res = sparse.gauss_seidel(&b, &x0, eps, max_iter);
    let seidel_sparse_elapsed = seidel_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(seidel_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_start = Instant::now();
    let gauss_res = mat.gaussian(&b);
    let gauss_elapsed = gauss_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gauss_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_res = match gauss_res {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };

    let gauss_sparse_start = Instant::now();
    let gauss_sparse_res = sparse.gaussian(&b);
    let gauss_sparse_elapsed = gauss_sparse_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(gauss_sparse_start.elapsed().subsec_nanos()) / 1_000_000.0;

    let gauss_sparse_res = match gauss_sparse_res {
        Ok(values) => values,
        Err(e) => {
            eprintln!("{}", e);
            process::exit(0);
        }
    };

    println!("\nMC:");
    let mc_start = Instant::now();
    let mc_res = monte_carlo::simulate_park_walk(&config, 30_000);
    let mc_elapsed = mc_start.elapsed().as_secs_f64() * 1000.0
        + f64::from(mc_start.elapsed().subsec_nanos()) / 1_000_000.0;
    println!("mc: {} in {:.6}ms", mc_res, mc_elapsed);

    println!("\nMAT:");
    println!(
        "gauss: {} in {:.6}ms",
        gauss_res[config.starting_pos], gauss_elapsed
    );
    println!(
        "gauss pp: {} in {:.6}ms",
        gpp_result[config.starting_pos], gpp_elapsed
    );
    println!(
        "gauss seidel: {} in {:.6}ms",
        seidel_res[config.starting_pos], seidel_elapsed
    );
    println!(
        "jacobi: {} in {:.6}ms",
        jacobi_res[config.starting_pos], jacobi_elapsed
    );
    println!("\nSPARSE:");
    println!(
        "gauss: {} in {:.6}ms",
        gauss_sparse_res[config.starting_pos], gauss_sparse_elapsed
    );
    println!(
        "gauss pp: {} in {:.6}ms",
        gpp_sparse_res[config.starting_pos], gpp_sparse_elapsed
    );
    println!(
        "gauss seidel: {} in {:.6}ms",
        seidel_sparse_res[config.starting_pos], seidel_sparse_elapsed
    );
    println!(
        "jacobi: {} in {:.6}ms",
        jacobi_sparse_res[config.starting_pos], jacobi_sparse_elapsed
    );
}
