use super::{ILP, Vector, ILPError, IntData};
use std::time::Instant;
use std::cmp::max;
use std::f64;

type Map<K,V> = hashbrown::HashMap<K,V>;
type LookupTable = Map<Vector, (Vector, i32)>;

/*
    based on https://arxiv.org/abs/1803.04744
*/

pub fn solve(ilp:&ILP) -> Result<Vector, ILPError> {
    println!("Solving ILP with the Discrepancy Algorithm...");
    let start = Instant::now();

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
    let b_bound = 4*H;

    println!(" -> H = {} >= herdisc(A)", H);
    println!(" -> Iterations: K = {}", K);

    let mut solutions     = LookupTable::with_capacity(ilp.A.num_cols());
    let mut new_solutions = LookupTable::with_capacity(ilp.A.num_cols());
    let mut x_bound = 1.0; // (6/5)^i
    let mut sb:Vector;  // scaled b (by 2^{i-K})

    // i=0 (trivial solutions)
    for (i, (column, &cost)) in ilp.A.iter().zip(ilp.c.iter()).enumerate() {
         solutions.insert(column.clone(), (Vector::unit(n, i), cost));
    }

    // i={1,...,K}
    for i in 1..K+1 {
        x_bound *= 1.2;
        let x_ibound = x_bound as IntData;
        sb = compute_sb(&ilp.b, K, i);
        
        // generate new solutions
        for (j, (b1, (x1,c1))) in solutions.iter().enumerate() {
            for (b2, (x2,c2)) in solutions.iter().skip(j+1) {
                let b = b1.add(b2);
                let x = x1.add(x2);
                let c = c1+c2;

                if x.one_norm() > x_ibound || !sb.max_distance(&b, b_bound) {
                    continue;
                }

                let insert = match new_solutions.get(&b) {
                    Some(&(_,cost)) => { cost < c },
                    None => true
                };

                if insert {
                    new_solutions.insert(b, (x, c));
                }
            }
        }

        for (b,x) in new_solutions.drain() {
            solutions.insert(b,x);
        }

        if i%50 == 0 {
            println!(" -> Iteration {}, size: {}, t: {:?}", i, solutions.len(), start.elapsed());
        }
    }

    println!(" -> {} iterations completed. t: {:?}", K, start.elapsed());

    match solutions.get(&ilp.b) {
        Some((x,_)) => Ok(x.clone()),
        None => Err(ILPError::NoSolution)
    }
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

fn compute_sb(b:&Vector, k:i32, i:i32) -> Vector {
    assert!(k>=i);
    let s = 0.5f64.powi(k-i);
    let mut v = Vector::new(b.len());

    for &bv in b.iter() {
        let x = bv as f64 * s;
        v.data.push(x.round() as IntData);
    }

    v
}