pub mod matrix;
pub mod monte_carlo;

fn main() {
    let mat = matrix::Matrix::from_vecs(vec![
        vec![1.2, 2.6, -0.1, 1.5],
        vec![4.5, 9.8, -0.4, 5.7],
        vec![0.1, -0.1, -0.3, -3.5],
        vec![4.5, -5.2, 4.2, -3.4],
    ]);
    let b = vec![13.15, 49.84, -14.08, -46.51];

    match mat.gaussian(&b) {
        Ok(solutions) => {
            println!("Gaussian solutions: {:.6?}", solutions);
        }
        Err(e) => {
            eprintln!("Error: {}", e);
        }
    }

    match mat.gaussian_partial_pivot(&b) {
        Ok(solutions) => {
            println!("Gaussian partial pivoting solutions: {:.6?}", solutions);
        }
        Err(e) => {
            eprintln!("Error: {}", e);
        }
    }
}
