use fnv::{FnvHashMap};
use std::ops::Range;
use super::Vector;
use std::slice::Iter;

type Edge = (usize,usize,i32,u16);

pub struct Node {
    idx: usize,
    edges: Vec<Edge>
}

pub struct VectorDiGraph {
    nodes: Vec<Node>,
    map: FnvHashMap<Vector, usize> // maps Vector -> idx
}

impl VectorDiGraph {
    pub fn new() -> Self {
        VectorDiGraph {
            nodes: Vec::new(),
            map: FnvHashMap::default()
        }
    }

    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    fn next_idx(&self) -> usize {
        self.size()
    }

    pub fn get_idx_by_vec(&self, v:&Vector) -> Option<usize> {
        match self.map.get(v) {
            Some(&idx) => Some(idx),
            None       => None
        }
    }

    pub fn add_node(&mut self, v:Vector) -> usize {
        let node = Node {
            idx: self.next_idx(),
            edges: Vec::new()
        };
        let node_idx = node.idx;
        self.nodes.push(node);
        self.map.insert(v, node_idx);

        node_idx
    }

    pub fn add_edge(&mut self, from: usize, to: usize, cost: i32, idx: u16) {
        let edge = (from, to, cost, idx);
        self.nodes[from].edges.push(edge);
    }

    pub fn iter_nodes(&self) -> Range<usize> {
        1..self.nodes.len()
    }

    pub fn iter_edges(&self, from_idx:usize) -> Iter<Edge> {
        self.nodes[from_idx].edges.iter()
    }
}
