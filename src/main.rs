pub mod ilp;
use self::ilp::*;

fn main() {
    println!("Hello, world!");

    let mat = Matrix::from_slice(3, 3, &[1,0,0,0,2,0,0,0,1]);
    let b = Vector::from_slice(&[5, 6, 5]);
    let c = Vector::from_slice(&[1, 2, 3]);

    let ilp = ILP::new(mat, b, c);

    match steinitz::solve(&ilp) {
        Ok(x) => println!("Found a solution! x={:?}", x),
        Err(ILPError::NoSolution) => println!("The ILP has no solution."),
        Err(ILPError::Unbounded)  => println!("The ILP is unbounded.")
    }
}

