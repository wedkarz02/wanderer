use std::io::{BufRead, BufReader};
use std::process::Command;
use std::time::Instant;
use std::{env, fs, process};

use base::*;
use comparisons::{compare_vecs, incremental_verify_mc};
use matrix::*;
use sparse::Sparse;

pub mod base;
pub mod comparisons;
pub mod matrix;
pub mod monte_carlo;
pub mod sparse;

#[derive(Debug)]
pub struct Sets(Vec<Vec<Vec<usize>>>);

#[derive(Debug, Clone)]
pub struct Intersection {
    pub id: usize,
    pub start: bool,
    pub well: bool,
    pub exit: bool,
}

impl Intersection {
    pub fn new(id: usize, start: bool, well: bool, exit: bool) -> Self {
        Self {
            id,
            start,
            well,
            exit,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Alley {
    pub a: Intersection,
    pub b: Intersection,
    pub length: usize,
}

impl Alley {
    pub fn new(a: Intersection, b: Intersection, length: usize) -> Self {
        Self { a, b, length }
    }

    pub fn get_propability(&self) -> f64 {
        let matrix_length = self.length + 2;
        let starting_pos = matrix_length - 2;
        1f64 - (starting_pos as f64 / (matrix_length - 1) as f64)
    }
}

#[derive(Debug, Clone)]
pub struct Config {
    pub inters: Vec<Intersection>,
    pub alleys: Vec<Alley>,
    pub starting_pos: usize,
}

impl Config {
    pub fn build(sets: Sets) -> Self {
        let inters_count = sets.0[0][0][0];
        let alleys_count = sets.0[0][0][1];

        let mut inters = Vec::new();
        for i in 1..inters_count + 1 {
            inters.push(Intersection::new(i, false, false, false));
        }

        for i in 1..sets.0[1][0].len() {
            for inter in &mut inters {
                if inter.id == sets.0[1][0][i] {
                    inter.well = true;
                }
            }
        }

        for i in 1..sets.0[1][1].len() {
            for inter in &mut inters {
                if inter.id == sets.0[1][1][i] {
                    inter.exit = true;
                }
            }
        }

        let mut starting_pos = 0;
        let start = sets.0[1][2][1];
        for inter in &mut inters {
            if inter.id == start {
                inter.start = true;
                starting_pos = inter.id - 1;
            }
        }

        let mut alleys = Vec::new();
        for i in 1..alleys_count + 1 {
            let a_id = sets.0[0][i][0];
            let b_id = sets.0[0][i][1];
            let mut a = Intersection::new(0, false, false, false);
            let mut b = Intersection::new(0, false, false, false);
            let length = sets.0[0][i][2];
            for inter in &inters {
                if inter.id == a_id {
                    a = inter.clone();
                }
                if inter.id == b_id {
                    b = inter.clone();
                }
            }
            alleys.push(Alley::new(a, b, length));
        }

        Self {
            inters,
            alleys,
            starting_pos,
        }
    }
}

pub fn parse_config(file_name: &'static str) -> Sets {
    let file = fs::File::open(file_name).expect("failed to open the file");
    let reader = BufReader::new(file);

    let mut data_sets: Vec<Vec<Vec<usize>>> = Vec::new();
    let mut current_set: Vec<Vec<usize>> = Vec::new();

    for line in reader.lines() {
        let line = line.expect("failed to read the line");
        if line.trim().is_empty() {
            data_sets.push(current_set);
            current_set = Vec::new();
        } else {
            let elements: Vec<usize> = line
                .split_whitespace()
                .map(|x| x.parse().expect("failed to parse the number"))
                .collect();

            current_set.push(elements);
        }
    }

    if !current_set.is_empty() {
        data_sets.push(current_set);
    }

    Sets(data_sets)
}

pub fn gen_config(inter_count: usize, alley_count: usize) -> Result<(), std::io::Error> {
    let py_output = Command::new("python3")
        .arg("scripts/gen_config.py")
        .arg(inter_count.to_string())
        .arg(alley_count.to_string())
        .output()
        .expect("failed to execute python process");

    if py_output.status.success() {
        let result = String::from_utf8_lossy(&py_output.stdout);
        fs::write("tmp.config", result.to_string())
    } else {
        let error = String::from_utf8_lossy(&py_output.stderr);
        eprintln!("{}", error);
        Ok(())
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("Nothing to do");
        process::exit(0);
    }

    match args[1].as_str() {
        "compare-config" => {
            let sets = parse_config(".config");
            let config = Config::build(sets);
            if let Err(e) = comparisons::incremental_compare_config(1, Some(&config)) {
                eprintln!("{}", e);
                process::exit(0);
            }
        }
        "compare-inc" => {
            if let Err(e) = comparisons::incremental_compare_config(100, None) {
                eprintln!("{}", e);
                process::exit(0);
            }

            let py_output = Command::new("python3")
                .arg("scripts/plot_config_cmps.py")
                .output()
                .expect("failed to execute python process");

            if py_output.status.success() {
                let result = String::from_utf8_lossy(&py_output.stdout);
                println!("{}", result);
            } else {
                let error = String::from_utf8_lossy(&py_output.stderr);
                eprintln!("{}", error);
            }
        }
        "compare-default" => {
            comparisons::incremental_compare_default();

            let py_output = Command::new("python3")
                .arg("scripts/plot_default_cmps.py")
                .output()
                .expect("failed to execute python process");

            if py_output.status.success() {
                let result = String::from_utf8_lossy(&py_output.stdout);
                println!("{}", result);
            } else {
                let error = String::from_utf8_lossy(&py_output.stderr);
                eprintln!("{}", error);
            }
        }
        "gen-config" => {
            if let Err(e) = gen_config(1000, 1500) {
                eprintln!("{}", e);
                process::exit(0);
            }
        }
        "check" => {
            let sets = parse_config("default.config");
            let config = Config::build(sets);

            let (mat, b) = Matrix::from_config(&config);
            let (sparse, _) = Sparse::from_config(&config);
            let x0 = vec![0f64; b.len()];
            let eps = 1e-16;
            let max_iter = 1_000_000;

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

            if compare_vecs(&gpp_sub, &b, 1e-16) {
                println!("Gauss (partial pivot): Success")
            } else {
                println!("Gauss (partial pivot): Failure")
            }
            if compare_vecs(&gpp_sparse_sub, &b, 1e-16) {
                println!("Sparse Gauss (partial pivot): Success")
            } else {
                println!("Sparse Gauss (partial pivot): Failure")
            }
            if compare_vecs(&jacobi_sub, &b, 1e-16) {
                println!("Jacobi: Success")
            } else {
                println!("Jacobi: Failure")
            }
            if compare_vecs(&jacobi_sparse_sub, &b, 1e-16) {
                println!("Sparse Jacobi: Success")
            } else {
                println!("Sparse Jacobi: Failure")
            }
            if compare_vecs(&seidel_sub, &b, 1e-16) {
                println!("Gauss-Seidel: Success")
            } else {
                println!("Gauss-Seidel: Failure")
            }
            if compare_vecs(&gauss_sub, &b, 1e-16) {
                println!("Gauss: Success")
            } else {
                println!("Gauss: Failure")
            }
            if compare_vecs(&gauss_sparse_sub, &b, 1e-16) {
                println!("Sparse Gauss: Success")
            } else {
                println!("Sparse Gauss: Failure")
            }
        }
        "time-all" => {
            let sets = parse_config("tmp.config");
            let config = Config::build(sets);
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
        "from-cfg" => {
            // NOTE:
            //   - Jacobi    -> OK
            //   - Gauss PP  -> OK
            //   - Gauss     -> OK
            //   - Seidel    -> OK
            //
            //   - Sparse Jacobi    -> OK
            //   - Sparse Gauss PP  -> Results OK but way to slow
            //   - Sparse Gauss     -> Results OK but way to slow
            //   - Sparse Seidel    -> Results OK but way to slow

            let sets = parse_config("tmp.config");
            let config = Config::build(sets);
            let (mat, _) = Matrix::from_config(&config);
            let (sparse, b) = Sparse::from_config(&config);

            let mc_res = monte_carlo::simulate_park_walk(&config, 10_000);
            println!("mc: {}", mc_res);

            let mat_res = mat.gaussian(&b).unwrap();
            println!("gauss: {:?}", mat_res[config.starting_pos]);

            let sp_res = sparse.gaussian(&b).unwrap();
            println!("sp gauss: {:?}", sp_res[config.starting_pos]);
        }
        "verify-mc" => {
            if let Err(e) = incremental_verify_mc(100, None) {
                eprintln!("{}", e);
                process::exit(0);
            }

            let py_output = Command::new("python3")
                .arg("scripts/plot_verify_mc.py")
                .output()
                .expect("failed to execute python process");

            if py_output.status.success() {
                let result = String::from_utf8_lossy(&py_output.stdout);
                println!("{}", result);
            } else {
                let error = String::from_utf8_lossy(&py_output.stderr);
                eprintln!("{}", error);
            }
        }
        _ => eprintln!("Unrecognised optional argument"),
    }
}
