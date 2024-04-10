pub mod matrix;

fn main() {
    let mut mat = matrix::Matrix::from_size(5, 10);

    for i in 0..mat.rows.len() {
        mat.rows[i][i] = 1.0;
    }

    println!("{}", mat);
}
