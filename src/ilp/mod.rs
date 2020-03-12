use std::cmp::max;
use std::slice::Iter;

pub mod steinitz;
mod graph;

pub type IntData = i32;

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
pub struct Vector {
    data: Vec<IntData>
}

pub struct Matrix {
    columns: Vec<Vector>,
    size: (usize, usize) // rows, columns or (m,n)
}

#[allow(non_snake_case)]
pub struct ILP {
    A: Matrix,
    b: Vector,
    c: Vector,
    delta: IntData,
    col_delta: IntData
}

pub enum ILPError {
    NoSolution,
    Unbounded
}

impl ILP {
    pub fn new(mat:Matrix, b:Vector, c:Vector) -> Self {
        assert!(b.len() == mat.size.0);
        assert!(c.len() == mat.size.1);
        assert!(mat.size.0 > 0 && mat.size.1 > 0);
    
        let mut delta = b.inf_norm();
    
        for v in mat.columns.iter() {
            delta = max(delta, v.inf_norm());
        }
    
        ILP {
            A: mat,
            b: b,
            c: c,
            delta: delta,
            col_delta: delta
        }
    }
}

impl Vector {
    pub fn new(size:usize) -> Self {
        Vector {
            data: Vec::with_capacity(size)
        }
    }

    pub fn zero(size:usize) -> Self {
        Vector {
            data: vec![0; size]
        }
    }

    pub fn from_slice(data:&[IntData]) -> Self {
        Vector {
            data: data.iter().cloned().collect()
        }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn iter(&self) -> Iter<IntData> {
        self.data.iter()
    }

    pub fn add(&self, other:&Vector) -> Vector {
        assert_eq!(self.len(), other.len());
        let mut v = Vec::with_capacity(self.len());

        for (x1,x2) in self.iter().zip(other.iter()) {
            v.push(x1 + x2);
        }

        Vector {
            data: v
        }
    }

    pub fn dot(&self, other: &Vector) -> IntData {
        assert_eq!(self.len(), other.len());
        let mut sum = 0;

        for (x1,x2) in self.iter().zip(other.iter()) {
            sum += x1*x2;
        }
    
        sum
    }

    pub fn norm2(&self) -> IntData {
        let mut sum = 0;

        for x in self.iter() {
            sum += x*x;
        }

        sum
    }

    pub fn norm(&self) -> f32 {
        let x = self.norm2() as f32;
        x.sqrt()
    }

    pub fn inf_norm(&self) -> IntData {
        let mut max = self.data[0];

        for &x in self.iter().skip(1) {
            if x > max {
                max = x;
            }
        }

        max
    }

    pub fn one_norm(&self) -> IntData {
        let mut sum = 0;

        for x in self.iter() {
            sum += x.abs();
        }

        sum
    }

    pub fn as_f32_vec(&self) -> Vec<f32> {
        let mut v = Vec::with_capacity(self.data.len());

        for &x in self.iter() {
            v.push(x as f32);
        }

        v
    }
}

impl Matrix {
    pub fn from_slice(rows:usize, columns:usize, data:&[IntData]) -> Matrix {
        assert_eq!(data.len(), rows*columns);
        let mut cols = Vec::with_capacity(columns);

        for c_idx in 0..columns {
            let start_idx = c_idx*rows;
            let end_idx = start_idx + rows;
            let c = Vector::from_slice(&data[start_idx..end_idx]);
            cols.push(c);
        }

        Matrix {
            columns: cols,
            size: (rows, columns)
        }
    }
}
