use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::{env, fs, process};

use base::*;
use comparisons::compare_vecs;
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
    pub trashcan: bool,
    pub start: bool,
    pub well: bool,
    pub exit: bool,
}

impl Intersection {
    pub fn new(id: usize, trashcan: bool, start: bool, well: bool, exit: bool) -> Self {
        Self {
            id,
            trashcan,
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
    pub deadends: Vec<usize>,
    pub starting_pos: usize,
}

impl Config {
    // This is not at all resistant to incorrect config file and will panic
    // left and right. It's not that important right now but should be fixed later.
    pub fn build(sets: Sets) -> Self {
        let inters_count = sets.0[0][0][0];
        let alleys_count = sets.0[0][0][1];

        let mut inters = Vec::new();
        for i in 1..inters_count + 1 {
            inters.push(Intersection::new(i, false, false, false, false));
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

        for i in 1..sets.0[1][3].len() {
            for inter in &mut inters {
                if inter.id == sets.0[1][3][i] {
                    inter.trashcan = true;
                }
            }
        }

        let mut alleys = Vec::new();
        for i in 1..alleys_count + 1 {
            let a_id = sets.0[0][i][0];
            let b_id = sets.0[0][i][1];
            let mut a = Intersection::new(0, false, false, false, false);
            let mut b = Intersection::new(0, false, false, false, false);
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

        let mut dead_ends: Vec<usize> = Vec::new();
        let mut not_dead_ends: Vec<usize> = Vec::new();

        for alley in &alleys {
            if !not_dead_ends.contains(&alley.a.id) {
                if dead_ends.contains(&alley.a.id) {
                    not_dead_ends.push(alley.a.id);
                    for i in 0..dead_ends.len() {
                        if dead_ends[i] == alley.a.id {
                            dead_ends.remove(i);
                            break;
                        }
                    }
                } else {
                    dead_ends.push(alley.a.id);
                }
            }
            if !not_dead_ends.contains(&alley.b.id) {
                if dead_ends.contains(&alley.b.id) {
                    not_dead_ends.push(alley.b.id);
                    for i in 0..dead_ends.len() {
                        if dead_ends[i] == alley.b.id {
                            dead_ends.remove(i);
                            break;
                        }
                    }
                } else {
                    dead_ends.push(alley.b.id);
                }
            }
        }

        Self {
            inters,
            alleys,
            deadends: dead_ends,
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
            if let Err(e) = gen_config(4, 4) {
                eprintln!("{}", e);
                process::exit(0);
            }
        }
        "verify" => {
            let sets = parse_config("default.config");
            let config = Config::build(sets);

            let mut map: HashMap<usize, f64> = HashMap::new();

            for alley in &config.alleys {
                map.insert(alley.length, alley.get_propability());
            }

            let max_iter = 10_000;
            let eps = 1e-2;

            for entry in map {
                let n = entry.0 + 2;
                let start = n - 2;
                let mc_result = monte_carlo::simulate_walk(n, start, max_iter);

                println!("built: {}\nmc: {}", entry.1, mc_result);
                if (entry.1 - mc_result).abs() > eps {
                    println!("Failed.\n");
                } else {
                    println!("Succeeded.\n");
                }
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
        "sim-walk" => {
            let sets = parse_config("default.config");
            let config = Config::build(sets);
            let (mat, b) = Matrix::from_config(&config);
            println!("{}{:?}", mat, b);

            let mc_res = monte_carlo::simulate_park_walk(&config, 10_000);
            println!("mc: {}", mc_res);

            let mat_res = mat.gaussian_partial_pivot(&b).unwrap()[config.starting_pos];
            println!("gs: {}", mat_res);
        }
        "asdf" => {
            let mat = Matrix::from_vecs(vec![
                vec![1.0, 0.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0, 0.0],
                vec![-6.0 / 14.0, -4.0 / 14.0, 1.0, -4.0 / 14.0],
                vec![0.0, 0.0, 0.0, 1.0],
            ]);
            let b = vec![0.0, 1.0, 0.0, 1.0];
            // let sets = parse_config("default.config");
            // let config = Config::build(sets);
            // let (mat, b) = Matrix::from_config(&config);
            println!("{}{:?}", mat, b);
            let mat_res = mat.gaussian_partial_pivot(&b).unwrap();
            println!("{:?}", mat_res);

            let bnew = mat.multiply_by_vec(&mat_res);
            println!("{:?}", bnew);
        }
        "from-cfg" => {
            let sets = parse_config("default.config");
            let config = Config::build(sets);
            let (mat, b) = Matrix::from_config(&config);

            println!("{}", mat);
            println!("{:?}", b);

            let mat_res = mat.gaussian_partial_pivot(&b).unwrap();
            println!("gs: {:?}", mat_res);

            let mc_res = monte_carlo::simulate_park_walk(&config, 10_000);
            println!("mc: {}", mc_res);

            let bnew = mat.multiply_by_vec(&mat_res);
            println!("{:?}", bnew);
        }
        _ => eprintln!("Unrecognised optional argument"),
    }
}
