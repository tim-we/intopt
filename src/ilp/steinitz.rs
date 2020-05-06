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
    println!("Solving ILP with the Eisenbrand & Weismantel algorithm...");
    let start = Instant::now();

    // constants
    let r = 1.0 / ilp.b.norm2() as f32;
    let (rows, columns) = ilp.A.size; // (m,n)
    let b_float = ilp.b.as_f32_vec();

    // graph
    let mut graph = VectorDiGraph::with_capacity(16384, columns);

    // construction surface
    let mut surface:Vec<(Vector, NodeIdx)> = Vec::with_capacity(16384);
    let mut new_surface:Vec<(Vector, NodeIdx)> = Vec::with_capacity(16384);
    let mut max_surface_size = 1;

    // add origin
    {
        let zero = Vector::zero(rows);
        graph.add_node(zero.clone(), 0, 0, 0);
        surface.push((zero, 0));
    }

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
        new_surface.reserve(max_new_nodes);

        // grow graph
        depth = depth+1;
        bound = compute_bound(ilp, depth);
        for (x, node_idx) in surface.drain(0..surface.len()) {
            let from = graph.get(node_idx).clone();

            // iterate over matrix columns
            for (i, (v,&c)) in ilp.A.iter().zip(ilp.c.iter()).enumerate() {
                // potentially new point
                let xp = x.add(v);
                let s = clamp(xp.dot(&ilp.b) as f32 * r, 0.0, 1.0);

                // ||xp - d*b|| <= bound
                if is_in_bounds(&xp, &b_float, s, bound) {
                    let cost = c as Cost;
                    let to_cost = from.cost + cost;

                    let to_idx = match graph.get_node_by_vec_mut(&xp) {
                        Some(node) => {
                            // this vector was already in the graph

                            // bellman-ford update 
                            if to_cost > node.cost {
                                node.predecessor = from.idx;
                                node.cost = to_cost;
                                node.via = i as ColumnIdx;
                            }

                            node.idx
                        },
                        None => {
                            // add new node
                            let idx = graph.add_node(xp.clone(), from.idx, to_cost, i as ColumnIdx);
                            new_surface.push((xp, idx));
                            idx
                        }
                    };

                    graph.add_edge(from.idx, to_idx, i as ColumnIdx);
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

    let b_node = match graph.get_node_by_vec(&ilp.b) {
        Some(node) => node.clone(),
        None => return Err(ILPError::NoSolution)
    };

    println!(" -> Continue Bellman-Ford Algorithm to find longest path...");
    let mut iterations = 0;
    // scan up to |V| - 2 times
    for _ in 2..graph.size() {
        let mut changed = false;
        iterations += 1;

        for node_idx in graph.iter_nodes() {
            let node = graph.get(node_idx).clone();
            for &(to, column) in node.edges.iter() {
                let to_cost = node.cost + ilp.c.data[column];
                let to_node = graph.get_mut(to);

                if to_cost > to_node.cost {
                    to_node.predecessor = node.idx;
                    to_node.cost = to_cost;
                    to_node.via = column;

                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }

    println!(" -> {} Bellman-Ford iterations, t={:?}", iterations, start.elapsed());
    println!(" -> Longest path cost: {}", b_node.cost);

    // create solution vector
    println!(" -> Creating solution vector... t={:?}", start.elapsed());

    let mut x = Vector::zero(columns);
    let b_idx = b_node.idx;
    let mut node = graph.get_node_by_vec_mut(&ilp.b).unwrap();

    // start from b and go backwards to 0
    loop {
        let pre = node.predecessor;

        if pre == b_idx {
            return Err(ILPError::Unbounded);
        } else {
            // mark node as visited
            node.predecessor = b_idx;
        }

        x.data[node.via as usize] += 1;
        node = graph.get_mut(pre);

        if node.idx == 0 {
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
