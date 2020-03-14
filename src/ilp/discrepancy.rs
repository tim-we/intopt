use super::{ILP, Vector, ILPError};
//use super::graph::*;
//use std::time::Instant;

/*
    based on https://arxiv.org/abs/1803.04744
*/

pub fn solve(ilp:&ILP) -> Result<Vector, ILPError> {
    println!("Solving ILP with the Discrepancy Algorithm...");

    if ilp.A.has_duplicate_columns() {
        println!(" -> The matrix has duplicate columns!");
        return Err(ILPError::UnsupportedMatrix);
    }
    
    // constants
    #[allow(non_snake_case)]
    let H = ilp.A.herdisc().ceil() as i32;
    println!(" -> herdisc(A) <= {} = H", H);

    unimplemented!();
}
