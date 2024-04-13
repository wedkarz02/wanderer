pub mod matrix;
pub mod monte_carlo;

fn main() {
    let n = 6;
    let mat = matrix::Matrix::init_default_path(n);
    let mut b = vec![0f64; n];
    b[0] = 1f64;
    let guess = vec![0f64; n];
    let starting_position = 3;

    let probs = mat.jacobi(&b, &guess, 20);
    let result = probs[starting_position];
    println!("Jacobi result: {}", result);

    let mc_result = monte_carlo::simulate_walk(n, starting_position, 500);
    println!("Monte Carlo result: {}", mc_result);
}
