use base::*;
use comparisons::*;
use sparse::*;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use std::{env, fs, process};

pub mod base;
pub mod comparisons;
pub mod matrix;
pub mod monte_carlo;
pub mod sparse;

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

#[derive(Debug)]
pub struct Alley {
    pub a: Intersection,
    pub b: Intersection,
    pub length: usize,
}

impl Alley {
    pub fn new(a: Intersection, b: Intersection, length: usize) -> Self {
        Self { a, b, length }
    }
}

#[derive(Debug)]
pub struct Config {
    pub inters: Vec<Intersection>,
    pub alleys: Vec<Alley>,
}

impl Config {
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

        let start = sets.0[1][2][1];
        for inter in &mut inters {
            if inter.id == start {
                inter.start = true;
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

        Self { inters, alleys }
    }
}

pub struct Sets(Vec<Vec<Vec<usize>>>);

fn parse_config(file_name: &'static str) -> Sets {
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

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("Nothing to do");
        process::exit(0);
    }

    let n = 1000;
    let eps = 1e-6;
    let max_iter = 1_000_000;
    let mc_max_iter = 10_000;
    let starting_pos = n / 2;

    match args[1].as_str() {
        "read" => {
            let set = parse_config(".config");
            let config = Config::build(set);
            println!("{:#?}", config);
        }
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

            let sparse = Sparse::init_default_path(n);
            let start = Instant::now();
            let res = sparse.gaussian_partial_pivot(&b);
            let time = start.elapsed();

            println!("{:?}\n{:#?}", time, res);
        }
        "sparse-jacobi" => {
            let n = 1000;
            let eps = 1e-16;
            let max_iter = 100_000_000;
            let starting_pos = n / 2;
            let sparse = Sparse::init_default_path(n);
            let mut b = vec![0f64; n];
            b[0] = 1f64;
            let x0 = vec![0f64; n];

            let jacobi_sparse_start = Instant::now();
            let jacobi_sparse_result = sparse.jacobi(&b, &x0, eps, max_iter)[starting_pos];
            let jacobi_sparse_elapsed = jacobi_sparse_start.elapsed();
            println!(
                "Jacobi sparse: {} in {:?}",
                jacobi_sparse_result, jacobi_sparse_elapsed
            );
        }
        "sparse-seidel" => {
            let n = 100;
            let eps = 1e-6;
            let max_iter = 1_000_000;
            let starting_pos = n / 2;
            let sparse = Sparse::init_default_path(n);
            // let mat = Matrix::init_default_path(n);
            let mut b = vec![0f64; n];
            b[0] = 1f64;
            let x0 = vec![0f64; n];

            let seidel_sparse_start = Instant::now();
            let seidel_sparse_result = sparse.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
            // let seidel_sparse_result = mat.gauss_seidel(&b, &x0, eps, max_iter)[starting_pos];
            let seidel_sparse_elapsed = seidel_sparse_start.elapsed();
            println!(
                "seidel sparse: {:?} in {:?}",
                seidel_sparse_result, seidel_sparse_elapsed
            );
        }
        _ => eprintln!("Unrecognised optional argument"),
    }
}
