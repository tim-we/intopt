## ILP Solver

A Rust program that solves ILPs (integer linear programs) using the algorithms described in [1] and [2].
Unfortunately the algorithms (or this implementation) are not very fast, especially the second one.  
Keep the number of constraints low!

## Usage

Install Rust & Cargo: https://rustup.rs  
Compile: `cargo build --release`  
Run: `cargo run --release -- examples/3x3.ilp` or `target/release/intopt examples/3x3.ilp`

Output for `examples/3x3.ilp`:
```
Reading file examples/3x3.ilp...
Parsing file...

ILP details:
 -> constraints: 3
 -> variables:   3
    ["x1", "x2", "x3"]
 -> Δ    = 2
 -> ‖b‖∞ = 6
 -> Matrix A:
|   1    0    0 |
|   0    2    0 |
|   0    0    1 |

 -> b = [5, 6, 5]
 -> c = [1, 2, 3]

Solving ILP with the Eisenbrand & Weismantel algorithm...
 -> Constructing the graph.............................
 -> Graph constructed! t=1.172876ms
    #vertices: 1070, #edges: 2864
    depth: 29, max. surface size: 79
    radius: start=12 end=6.6206894
 -> Continue Bellman-Ford Algorithm to find longest path...
 -> 1 Bellman-Ford iterations, t=1.203325ms
 -> Longest path cost: 26
 -> Creating solution vector... t=1.210861ms
 -> Done! Time elapsed: 1.214739ms

Solution:
 x1 = 5
 x2 = 3
 x3 = 5
``` 

## Papers

[1] https://arxiv.org/abs/1707.00481v3  
[2] https://arxiv.org/abs/1803.04744v3
