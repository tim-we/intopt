use std::ops::Range;
use super::{Vector, Cost};

type Map<K,V> = hashbrown::HashMap<K,V>;
pub type NodeIdx = usize;
pub type ColumnIdx = usize;

/*  A node contains its outgoing edges, thus an edge only
    stores one end index and the column that was used.
 */
pub type Edge = (NodeIdx, ColumnIdx);

#[derive(Clone)]
pub struct Node {
    pub idx: NodeIdx,
    pub predecessor: NodeIdx,
    pub via: ColumnIdx,
    pub cost: Cost,
    pub edges: Vec<Edge>
}

pub struct VectorDiGraph {
    nodes: Vec<Node>,
    map: Map<Vector, NodeIdx>,
    edges_per_node: usize, // avoid re-allocation
    edges: usize
}

impl VectorDiGraph {
    pub fn with_capacity(node_capacity:usize, edges:usize) -> Self {
        VectorDiGraph {
            nodes: Vec::with_capacity(node_capacity),
            map: Map::with_capacity(node_capacity),
            edges_per_node: edges,
            edges: 0
        }
    }

    pub fn reserve(&mut self, additional: usize) {
        self.nodes.reserve(additional);
        self.map.reserve(additional);
    }

    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    fn next_idx(&self) -> NodeIdx {
        self.size() as NodeIdx
    }

    pub fn get(&self, idx:NodeIdx) -> &Node {
        &self.nodes[idx]
    }

    pub fn get_mut(&mut self, idx:NodeIdx) -> &mut Node {
        &mut self.nodes[idx]
    }

    pub fn get_node_by_vec(&self, v:&Vector) -> Option<&Node> {
        match self.map.get(v) {
            Some(&idx) => Some(&self.nodes[idx]),
            None       => None
        }
    }

    pub fn get_node_by_vec_mut(&mut self, v:&Vector) -> Option<&mut Node> {
        match self.map.get(v) {
            Some(&idx) => Some(&mut self.nodes[idx]),
            None       => None
        }
    }

    pub fn add_node(&mut self, v:Vector, pre:NodeIdx, cost:Cost, via:ColumnIdx) -> NodeIdx {
        let node = Node {
            idx: self.next_idx(),
            edges: Vec::with_capacity(self.edges_per_node),
            predecessor: pre,
            via: via,
            cost: cost
        };
        let node_idx = node.idx;
        self.nodes.push(node);
        self.map.insert(v, node_idx);

        node_idx
    }

    pub fn add_edge(&mut self, from: NodeIdx, to: NodeIdx, idx: ColumnIdx) {
        let edge = (to, idx);
        self.nodes[from].edges.push(edge);
        self.edges += 1;
    }

    pub fn iter_nodes(&self) -> Range<usize> {
        1..self.nodes.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges
    }
}
