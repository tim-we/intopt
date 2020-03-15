use super::{ILP, Vector, ILPError};
use std::time::Instant;
use std::cmp::max;

type Map<K,V> = hashbrown::HashMap<K,V>;
type LookupTable = Map<Vector, (Vector, i32)>;

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
    let (_,n) = ilp.A.size;
    #[allow(non_snake_case)]
    let H = ilp.A.herdisc_upper_bound().ceil() as i32;
    #[allow(non_snake_case)]
    let K = compute_K(ilp);

    println!(" -> herdisc(A) <= {} = H", H);
    println!(" -> Iterations: K = {}", K);

    let mut solutions     = LookupTable::with_capacity(ilp.A.num_cols());
    let mut new_solutions = LookupTable::with_capacity(ilp.A.num_cols());

    // i=0 (trivial solutions)
    for (i, (column, &cost)) in ilp.A.iter().zip(ilp.c.iter()).enumerate() {
         solutions.insert(column.clone(), (Vector::unit(n, i), cost));
    }

    // i={1,...,K}
    let mut _x_bound = 1.0;
    for _ in 0..K {
        _x_bound *= 1.2;
        
        // generate new solutions
        for (j, (b1, (x1,c1))) in solutions.iter().enumerate() {
            for (b2, (x2,c2)) in solutions.iter().skip(j+1) {
                let _b3 = b1.add(b2);
                let _x3 = x1.add(x2);
                let _c3 = c1+c2;
                unimplemented!();
            }
        }

        solutions.clear();

        // swap
        {
            let tmp = solutions;
            solutions = new_solutions;
            new_solutions = tmp;
        }

    }

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

