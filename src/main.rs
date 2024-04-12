use std::process;

pub mod matrix;

fn main() {
    let mut mat = matrix::Matrix::from_size(2, 3);
    let mut other = matrix::Matrix::from_size(3, 2);
    mat.rows[0][0] = 1f64;
    mat.rows[0][1] = 4f64;
    mat.rows[0][2] = 2f64;
    mat.rows[1][0] = 2f64;
    mat.rows[1][1] = 5f64;
    mat.rows[1][2] = 1f64;

    other.rows[0][0] = 2f64;
    other.rows[0][1] = 1f64;
    other.rows[1][0] = 1f64;
    other.rows[1][1] = 1f64;
    other.rows[2][0] = 3f64;
    other.rows[2][1] = 2f64;

    let res = mat.multiply(&other).unwrap_or_else(|err| {
        eprintln!("{}", err);
        process::exit(1);
    });

    println!("{}", res);
}
