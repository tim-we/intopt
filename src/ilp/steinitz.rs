use num_traits::Float;
use super::{ILP, Vector, ILPError};
use std::time::Instant;
use super::graph::*;

type Set<T> = hashbrown::HashSet<T>; //fnv::FnvHashSet<T>;

/* 
    based on https://arxiv.org/abs/1707.00481v3    

    vertices = {0}
    edges = {}
    surface = {0}
    new_surface = {}
    ib = 1.0/||b||
    while surface not empty
        for all columns v_i in A:
            for all x in surface
                x' = x + v_i
                h = clamp(<x',b>*ib, 0, 1) * b //closest point on line
                if ||x' - h|| <= bound 
                    if x' not in vertices
                        add x' to vertices
                        add x' to new_surface
                    add edge (x,x') to edges with weight c_i
        surface = new_surface
        new_surface = {}
    compute longest path
    PROFIT
*/
pub fn solve(ilp:&ILP) -> Result<Vector, ILPError> {
    println!("Solving ILP with the Steinitz Algorithm...");
    let start = Instant::now();

    // constants
    let r = 1.0 / ilp.b.norm() as f32;
    let (rows, columns) = ilp.A.size; // (m,n)
    let b_float = ilp.b.as_f32_vec();
    let bound = 1.5 * (rows as i32 * ilp.delta) as f32;
    let ts_size_bound = ((2.0*bound) as i32 + 1).pow(rows as u32);
    println!(" -> Using {} as the bound for the tube set.", bound);
    println!(" -> Tube set size bound: {}", ts_size_bound);

    // graph
    let mut graph = VectorDiGraph::with_capacity(16384, columns);

    // construction surface
    let mut surface:Set<Vector> = Set::with_capacity(16384);
    let mut new_surface:Set<Vector> = Set::with_capacity(16384);
    let mut max_surface_size = 1;

    // add origin
    {
        let zero = Vector::zero(rows);
        graph.add_node(zero.clone());
        surface.insert(zero);
    }

    // bellman-ford data (distance, predecessor, matrix column index)
    let mut bf_data = Vec::<(Cost, NodeIdx, ColumnIdx)>::new();
    bf_data.push((0,0,0));

    // construct graph
    println!(" -> Constructing the graph...");
    while !surface.is_empty() {
        for x in surface.drain() {
            let from_idx = graph.get_idx_by_vec(&x).unwrap();

            for (v,i) in ilp.A.columns.iter().zip(0..columns) {
                let xp = x.add(v);
                let s = clamp(xp.dot(&ilp.b) as f32 * r, 0.0, 1.0);

                // ||xp - d*b|| <= bound
                if is_in_bounds(&xp, &b_float, s, bound) {
                    let cost = ilp.c.data[i];
                    let to_distance = bf_data[from_idx].0 + cost;

                    let to_idx = match graph.get_idx_by_vec(&xp) {
                        Some(to_idx) => {
                            // this vector was already in the graph

                            // bellman-ford update 
                            if to_distance > bf_data[to_idx].0 {
                                bf_data[to_idx].0 = to_distance;
                                bf_data[to_idx].1 = from_idx;
                            }

                            to_idx
                        },
                        None => {
                            // add new node
                            new_surface.insert(xp.clone());
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

    println!(" -> Graph constructed! t={:?}", start.elapsed());
    println!("    size: {}, max. surface size: {}", graph.size(), max_surface_size);

    let b_idx = match graph.get_idx_by_vec(&ilp.b) {
        Some(idx) => idx,
        None => return Err(ILPError::NoSolution)
    };

    println!(" -> Continue Bellman-Ford Algorithm to find longest path...");
    let mut iterations = 0;
    // scan up to |V| - 1 times
    for _ in 1..graph.size() {
        let mut changed = false;
        iterations += 1;

        for node_idx in graph.iter_nodes() {
            for &(from, to, cost, _) in graph.iter_edges(node_idx) {
                let to_distance = bf_data[from].0 + cost;
                if to_distance > bf_data[to].0 {
                    bf_data[to].0 = to_distance;
                    bf_data[to].1 = from;
                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }

    println!(" -> {} Bellman-Ford iterations, t={:?}", iterations, start.elapsed());

    println!(" -> Longest path distance: {:?}", bf_data[b_idx].0);

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
    assert!(min <= max);

    T::min(T::max(min, x), max)
}

/// ||x - s*b||_{inf} <= bound
fn is_in_bounds(x:&Vector, b:&Vec<f32>, s:f32, bound:f32) -> bool {
    let mut max = 0.0;

    for i in 0..x.len() {
        let d = (x.data[i] as f32 - (s * b[i])).abs();
        
        if d > max {
            max = d;
        }
    }

    max <= bound
}
