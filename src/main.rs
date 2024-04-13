pub mod matrix;

fn main() {
    let mat = matrix::Matrix::from_vecs(vec![vec![2.0, 1.0], vec![5.0, 7.0]]);

    let b = vec![11.0, 13.0];
    let x = vec![1.0; b.len()];

    println!("{}\n{:?}\n\n{:?}", mat, b, x);
    println!("{:?}", mat.jacobi(&b, &x, 25));
}
