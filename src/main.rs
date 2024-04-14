use std::time::Instant;

pub mod matrix;
pub mod monte_carlo;

fn main() {
    let n = 100;
    let eps = 1e-3;
    let mat = matrix::Matrix::init_default_path(n);
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let x0 = vec![0.0; b.len()];

    let j_start = Instant::now();
    let jacbi_res = mat.jacobi(&b, &x0, eps, 10000);
    let j_elapsed = j_start.elapsed();

    let s_start = Instant::now();
    let seidel_res = mat.gauss_seidel(&b, &x0, eps, 10000);
    let s_elapsed = s_start.elapsed();

    let g_start = Instant::now();
    let g_res = mat.gaussian(&b).unwrap();
    let g_elapsed = g_start.elapsed();

    let gp_start = Instant::now();
    let gp_res = mat.gaussian_partial_pivot(&b).unwrap();
    let gp_elapsed = gp_start.elapsed();

    println!("jacobi elapsed:          {:?} {:.6?}", j_elapsed, jacbi_res);
    println!(
        "seidel elapsed:          {:?} {:.6?}",
        s_elapsed, seidel_res
    );
    println!(
        "gauss without pivot elapsed:  {:?} {:.6?}",
        g_elapsed, g_res
    );
    println!(
        "gauss partial pivot elapsed:  {:?} {:.6?}",
        gp_elapsed, gp_res
    );
}
