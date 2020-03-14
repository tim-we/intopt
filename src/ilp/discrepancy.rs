use super::{ILP, Vector, ILPError};
//use super::graph::*;
use std::time::Instant;
use std::cmp::max;

/*
    based on https://arxiv.org/abs/1803.04744
*/

pub fn solve(ilp:&ILP) -> Result<Vector, ILPError> {
    println!("Solving ILP with the Discrepancy Algorithm...");
    let _start = Instant::now();

    if ilp.A.has_duplicate_columns() {
        println!(" -> The matrix has duplicate columns!");
        return Err(ILPError::UnsupportedMatrix);
    }
    
    // constants
    #[allow(non_snake_case)]
    let H = ilp.A.herdisc_upper_bound().ceil() as i32;
    #[allow(non_snake_case)]
    let K = compute_K(ilp);

    println!(" -> herdisc(A) <= {} = H", H);
    println!(" -> Iterations: K = {}", K);

    unimplemented!();
}

#[allow(non_snake_case)]
fn compute_K(ilp:&ILP) -> i32 {
    let n = ilp.A.size.0 as f64;
    let m = ilp.A.size.0 as i32;

    let x1 = f64::ln((m*max(ilp.delta_A, ilp.delta_b)) as f64);
    let x2 = (2*m+1) as f64 * x1;
    let x3 = 2.0 * f64::ln(n);
    let x4 = f64::ln(1.2);

    f64::ceil((x3 + x2)/x4) as i32
}
