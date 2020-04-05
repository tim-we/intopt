use num_traits::Float;
use super::{ILP, Vector, ILPError, Cost};
use std::time::Instant;
use super::graph::*;
use std::io;
use std::io::Write;
use ignore_result::Ignore;

/* 
    based on https://arxiv.org/abs/1707.00481v3
*/

pub fn solve(ilp:&ILP) -> Result<Vector, ILPError> {
    println!("Solving ILP with the Steinitz Algorithm...");
    let start = Instant::now();

    if ilp.A.has_zero_columns() {
        println!(" -> Matrix contains zero columns!");
        return Err(ILPError::Unsupported);
    }

    // constants
    let r = 1.0 / ilp.b.norm();
    let (rows, columns) = ilp.A.size; // (m,n)
    let b_float = ilp.b.as_f32_vec();

    // graph
    let mut graph = VectorDiGraph::with_capacity(16384, columns);

    // construction surface
    let mut surface:Vec<Vector> = Vec::with_capacity(16384);
    let mut new_surface:Vec<Vector> = Vec::with_capacity(16384);
    let mut max_surface_size = 1;

    // add origin
    {
        let zero = Vector::zero(rows);
        graph.add_node(zero.clone());
        surface.push(zero);
    }

    // bellman-ford data (distance, predecessor, matrix column [index] that was used to get to this node)
    let mut bf_data = Vec::<(Cost, NodeIdx, ColumnIdx)>::with_capacity(16384);
    bf_data.push((0,0,0));

    // construct graph
    print!(" -> Constructing the graph");
    io::stdout().flush().ignore();

    let mut bound;
    let mut depth = 0;
    while !surface.is_empty() {
        print!(".");
        io::stdout().flush().ignore();
        
        // pre-allocate memory for new nodes
        let max_new_nodes = surface.len() * columns;
        graph.reserve(max_new_nodes);
        bf_data.reserve(max_new_nodes);
        new_surface.reserve(max_new_nodes);

        // grow graph
        depth = depth+1;
        bound = compute_bound(ilp, depth);
        for x in surface.drain(0..surface.len()) {
            let from_idx = graph.get_idx_by_vec(&x).unwrap();

            for (i, (v,&c)) in ilp.A.iter().zip(ilp.c.iter()).enumerate() {
                let xp = x.add(v);
                let s = clamp(xp.dot(&ilp.b) as f32 * r, 0.0, 1.0);

                // ||xp - d*b|| <= bound
                if is_in_bounds(&xp, &b_float, s, bound) {
                    let cost = c as Cost;
                    let to_distance = bf_data[from_idx].0 + cost;

                    let to_idx = match graph.get_idx_by_vec(&xp) {
                        Some(to_idx) => {
                            // this vector was already in the graph

                            // bellman-ford update 
                            if to_distance > bf_data[to_idx].0 {
                                bf_data[to_idx].0 = to_distance;
                                bf_data[to_idx].1 = from_idx;
                                bf_data[to_idx].2 = i as ColumnIdx;
                            }

                            to_idx
                        },
                        None => {
                            // add new node
                            new_surface.push(xp.clone());
                            bf_data.push((to_distance, from_idx, i as ColumnIdx));
                            graph.add_node(xp)
                        }
                    };

                    graph.add_edge(from_idx, to_idx, cost, i as ColumnIdx);
                }
            }
        }

        // swap buffers (keep capacity/avoid new allocation)
        {
            let tmp = surface;
            surface = new_surface;
            new_surface = tmp;
        }

        if surface.len() > max_surface_size {
            max_surface_size = surface.len();
        }
    }

    println!();
    println!(" -> Graph constructed! t={:?}", start.elapsed());
    println!("    #vertices: {}, #edges: {}", graph.size(), graph.num_edges());
    println!("    depth: {}, max. surface size: {}", depth, max_surface_size);
    println!("    radius: start={} end={}", compute_bound(ilp, 1), compute_bound(ilp, depth));

    let b_idx = match graph.get_idx_by_vec(&ilp.b) {
        Some(idx) => idx,
        None => return Err(ILPError::NoSolution)
    };

    println!(" -> Continue Bellman-Ford Algorithm to find longest path...");
    let mut iterations = 0;
    // scan up to |V| - 2 times
    for _ in 2..graph.size() {
        let mut changed = false;
        iterations += 1;

        for node_idx in graph.iter_nodes() {
            for &(from, to, cost, i) in graph.iter_edges(node_idx) {
                let to_distance = bf_data[from].0 + cost;
                if to_distance > bf_data[to].0 {
                    bf_data[to].0 = to_distance;
                    bf_data[to].1 = from;
                    bf_data[to].2 = i;

                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }

    println!(" -> {} Bellman-Ford iterations, t={:?}", iterations, start.elapsed());
    println!(" -> Longest path cost: {}", bf_data[b_idx].0);

    // create solution vector
    println!(" -> Creating solution vector... t={:?}", start.elapsed());

    let mut x = Vector::zero(columns);
    let mut node = b_idx;

    // start from b and go backwards to 0
    loop {
        let (_, predecessor, column) = bf_data[node];

        if predecessor == b_idx {
            return Err(ILPError::Unbounded);
        } else {
            // mark node as visited
            bf_data[node].1 = b_idx;
        }

        node = predecessor;
        x.data[column as usize] += 1;

        if node == 0 {
            break;
        }
    }

    println!(" -> Done! Time elapsed: {:?}", start.elapsed());

    Ok(x)
}

fn clamp<T: Float>(x:T, min: T, max: T) -> T {
    debug_assert!(min <= max);

    T::min(T::max(min, x), max)
}

fn compute_bound(ilp:&ILP, depth:i32) -> f32 {
    let (m,_) = ilp.A.size;
    let da = ilp.delta_A as f32;
    let db = ilp.delta_b as f32;
    let delta = f32::min(
        2.0 * da,
        da + (1.0/depth as f32) * db
    );
    delta * m as f32
}

/// ||x - s*b||_{inf} <= bound
fn is_in_bounds(v:&Vector, b:&Vec<f32>, s:f32, bound:f32) -> bool {
    debug_assert!(v.len() == b.len());

    for (&x,&b) in v.iter().zip(b.iter()) {
        let d = (x as f32 - (s * b)).abs();
        
        if d > bound {
            return false;
        }
    }

    true
}
