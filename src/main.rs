use std::io::{BufRead, BufReader};
use std::process::Command;
use std::{env, fs, process};

use base::*;
use comparisons::incremental_verify_mc;
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
    if args.len() < 2 {
        println!("Nothing to do");
        process::exit(0);
    }

    match args[1].as_str() {
        "compare-config" => {
            let sets = parse_config(".config");
            let config = Config::build(sets);
            if let Err(e) = comparisons::incremental_compare_config(1, 10, Some(&config)) {
                eprintln!("{}", e);
                process::exit(0);
            }
        }
        "compare-inc" => {
            if let Err(e) = comparisons::incremental_compare_config(100, 20, None) {
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
            if args.len() != 4 {
                println!("gen-config requires 2 'usize' arguments");
                process::exit(0);
            }

            let (inter_count, alley_count): (usize, usize) =
                match (args[2].parse(), args[3].parse()) {
                    (Ok(a), Ok(b)) => (a, b),
                    _ => {
                        eprintln!("parsing args failed");
                        process::exit(0);
                    }
                };

            if let Err(e) = gen_config(inter_count, alley_count) {
                eprintln!("{}", e);
                process::exit(0);
            }
        }
        "check-result" => {
            let sets = parse_config("default.config");
            let config = Config::build(sets);

            comparisons::check_results(&config);
        }
        "time-all" => {
            let sets = parse_config("tmp.config");
            let config = Config::build(sets);

            comparisons::time_all(&config);
        }
        "from-cfg" => {
            let sets = parse_config("tmp.config");
            let config = Config::build(sets);
            let (mat, _) = Matrix::from_config(&config);
            let (sparse, b) = Sparse::from_config(&config);

            let mc_res = monte_carlo::simulate_park_walk(&config, config.inters.len() * 10);
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
