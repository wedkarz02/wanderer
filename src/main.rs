pub mod matrix;
pub mod monte_carlo;

fn main() {
    // let n = 6;
    // let mat = matrix::Matrix::init_default_path(n);
    // let mut b = vec![0f64; n];
    // b[0] = 1f64;
    // let guess = vec![0f64; n];
    // let starting_position = 3;

    // let probs = mat.jacobi(&b, &guess, 20);
    // let result = probs[starting_position];
    // println!("Jacobi result: {}", result);

    // let mc_result = monte_carlo::simulate_walk(n, starting_position, 500);
    // println!("Monte Carlo result: {}", mc_result);

    // if let Ok(solution) = mat.gaussian(&b) {
    //     println!(
    //         "Gaussian elimination result: {:?}",
    //         solution[starting_position]
    //     );
    // }

    let mat = matrix::Matrix::from_vecs(vec![
        vec![1.2, 2.6, -0.1, 1.5],
        vec![4.5, 9.8, -0.4, 5.7],
        vec![0.1, -0.1, -0.3, -3.5],
        vec![4.5, -5.2, 4.2, -3.4],
    ]);
    let b = vec![13.15, 49.84, -14.08, -46.51];

    match mat.gaussian(&b) {
        Ok(solutions) => {
            println!("Gaussian solutions: {:?}", solutions);
        }
        Err(e) => {
            eprintln!("Error: {}", e);
        }
    }
}
