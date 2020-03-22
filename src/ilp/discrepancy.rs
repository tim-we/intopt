use super::{ILP, Vector, ILPError, IntData, Matrix, Cost};
use std::time::Instant;
use std::cmp::max;
use std::{f64, i32};

type Map<K,V> = hashbrown::HashMap<K,V>;
type LookupTable = Map<Vector, (Vector, Cost)>;

/*
    based on https://arxiv.org/abs/1803.04744
*/

pub fn solve(ilp:&ILP) -> Result<Vector, ILPError> {
    println!("Solving ILP with the Discrepancy Algorithm...");
    let start = Instant::now();

    if ilp.b.inf_norm() == 0 {
        println!(" -> b=0");
        return Err(ILPError::Unsupported);
    }

    let modified_ilp;

    let p = if ilp.A.has_duplicate_columns() {
            println!(" -> The matrix has duplicate columns!");
            let (new_ilp, removed) = ilp_with_unique_columns(ilp);
            println!(" -> Solving modified ILP (removed {} columns)", removed);
            new_ilp.print_details("    ");

            modified_ilp = new_ilp;
            &modified_ilp
        } else {
            ilp
        };
    
    // constants
    let (m,n) = p.A.size;
    #[allow(non_snake_case)]
    let H = p.A.herdisc_upper_bound().ceil() as i32;
    #[allow(non_snake_case)]
    let K = compute_K(ilp);
    let b_bound = 4*H;

    println!(" -> H = {} >= herdisc(A)", H);
    println!(" -> K = {}", K);

    let mut solutions = LookupTable::with_capacity(1024);
    
    // i=0 (trivial solutions)
    solutions.insert(Vector::zero(m), (Vector::zero(n), 0));
    for (i, (column, &cost)) in p.A.iter().zip(p.c.iter()).enumerate() {
        solutions.insert(column.clone(), (Vector::unit(n, i), cost));
    }

    // pre-compute main iteration (scaled b, max iterations, max x bound)
    let mut iterations = Vec::<(Vector, usize)>::with_capacity(max(p.delta_b as usize, 2));
    {
        let mut last = (compute_sb(&p.b, K, 1), 1); // i=1

        // i={1,...,K}
        for i in 1..K+1 {
            let sb = compute_sb(&p.b, K, i); // scaled b (by 2^{i-K})

            if sb != last.0 {
                iterations.push(last);
                last = (sb, 1);

                if i == K {
                    iterations.push(last.clone());
                }
            } else {
                last.1 += 1;
            }
        }

        assert_eq!(last.0, p.b);
        println!(" -> Iterations: {}", iterations.len());
    }

    let mut last_solutions = solutions.clone();
    let mut new_solutions  = LookupTable::with_capacity(512);
    
    println!(" -> Building lookup table...");
    for (sb, it_max) in iterations {
        println!("    > size: {}", solutions.len());

        for j in 0..it_max {
            // generate new solutions
            let iterator = if j==0 { solutions.iter() } else { last_solutions.iter() };
            for (k, (b1, (x1,c1))) in iterator.enumerate() {
                for (b2, (x2,c2))  in solutions.iter().skip(if j==0 {k+1} else {0}) {
                    let b = b1.add(b2);
                    let x = x1.add(x2);
                    let c = c1+c2;

                    if !sb.max_distance(&b, b_bound) {
                        continue;
                    }

                    let insert = match solutions.get(&b) {
                        Some(&(_,cost)) => { cost < c },
                        None => true
                    };

                    if insert {
                        new_solutions.insert(b, (x,c));
                    }
                }
            }

            // if there are no new solutions we can skip iterations j+1..it_max
            if new_solutions.is_empty() {
                continue;
            }

            // update lookup table
            for (b,x) in new_solutions.iter() {
                solutions.insert(b.clone(), x.clone());
            }

            // swap buffers
            {
                let tmp = last_solutions;
                last_solutions = new_solutions;
                new_solutions = tmp;
                new_solutions.clear();
            }
        }

        last_solutions.clear();
    }

    println!(" -> Done. Time elapsed: {:?}", start.elapsed());

    match solutions.get(&p.b) {
        Some((x,_)) => {
            println!(" -> Solution cost: {}", x.dot(&p.c));
            Ok(x.clone())
        },
        None => Err(ILPError::NoSolution)
    }
}

#[allow(non_snake_case)]
fn compute_K(ilp:&ILP) -> usize {
    let n = ilp.A.size.0 as f64;
    let m = ilp.A.size.0 as i32;

    let x1 = f64::ln((m*max(ilp.delta_A, ilp.delta_b)) as f64);
    let x2 = (2*m+1) as f64 * x1;
    let x3 = 2.0 * f64::ln(n);
    let x4 = f64::ln(1.2);

    f64::ceil((x3 + x2)/x4) as usize
}

fn compute_sb(b:&Vector, k:usize, i:usize) -> Vector {
    debug_assert!(k >= i);
    let s = 0.5f64.powi(k as i32 - i as i32);
    let mut v = Vector::new(b.len());

    for &bv in b.iter() {
        let x = bv as f64 * s;
        v.data.push(x.round() as IntData);
    }

    v
}

fn ilp_with_unique_columns(ilp:&ILP) -> (ILP, usize) {
    let mut mat = Matrix {
        columns: Vec::new(),
        size: (ilp.b.len(), 0)
    };

    let mut c = Vector {
        data: Vec::new()
    };

    let mut removed = Vec::new();
    for (i, col1) in ilp.A.iter().enumerate() {
        if removed.contains(&i) {
            continue;
        }

        let mut best = (col1, ilp.c.data[i]);
        for (j, col2) in ilp.A.iter().enumerate().skip(i+1) {
            if col1 == col2 {
                let cost = ilp.c.data[j];
                // keep column with highest cost/weight
                if cost > best.1 {
                    best = (col2, cost);
                    removed.push(i);
                } else {
                    removed.push(j);
                }
            }
        }
        
        mat.columns.push(best.0.clone());
        mat.size.1 += 1;
        c.data.push(best.1);
    }

    (ILP::new(mat, ilp.b.clone(), c), removed.len())
}
