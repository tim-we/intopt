use fnv::{FnvHashMap};
use std::ops::Range;
use super::Vector;
use std::slice::Iter;

pub type NodeIdx = usize;
pub type Cost = i32;
pub type ColumnIdx = u8;

pub type Edge = (NodeIdx, NodeIdx, Cost, ColumnIdx);

pub struct Node {
    idx: NodeIdx,
    edges: Vec<Edge>
}

pub struct VectorDiGraph {
    nodes: Vec<Node>,
    map: FnvHashMap<Vector, NodeIdx>
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

    fn next_idx(&self) -> NodeIdx {
        self.size() as NodeIdx
    }

    pub fn get_idx_by_vec(&self, v:&Vector) -> Option<NodeIdx> {
        match self.map.get(v) {
            Some(&idx) => Some(idx),
            None       => None
        }
    }

    pub fn add_node(&mut self, v:Vector) -> NodeIdx {
        let node = Node {
            idx: self.next_idx(),
            edges: Vec::new()
        };
        let node_idx = node.idx;
        self.nodes.push(node);
        self.map.insert(v, node_idx);

        node_idx
    }

    pub fn add_edge(&mut self, from: NodeIdx, to: NodeIdx, cost: Cost, idx: ColumnIdx) {
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
