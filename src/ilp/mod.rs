use std::cmp::max;

pub mod steinitz;
mod graph;

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
pub struct Vector {
    data: Vec<i32>
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
    delta: i32
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
    
        let mut delta:i32 = b.inf_norm();
    
        for v in mat.columns.iter() {
            delta = max(delta, v.inf_norm());
        }
    
        ILP {
            A: mat,
            b: b,
            c: c,
            delta: delta
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

    pub fn from_slice(data:&[i32]) -> Self {
        Vector {
            data: data.iter().cloned().collect()
        }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn add(&self, other:&Vector) -> Vector {
        assert_eq!(self.len(), other.len());
        let mut v = Vec::with_capacity(self.len());

        for i in 0..self.len() {
            v.push(self.data[i] + other.data[i]);
        }

        Vector {
            data: v
        }
    }

    pub fn dot(&self, other: &Vector) -> i32 {
        assert_eq!(self.len(), other.len());
    
        let mut sum = 0;
        for i in 0..self.len() {
            sum += self.data[i] * other.data[i];
        }
    
        sum
    }

    pub fn norm2(&self) -> i64 {
        let mut sum = 0i64;

        for i in 0..self.len() {
            sum += self.data[i] as i64 * self.data[i] as i64;
        }

        sum
    }

    pub fn norm(&self) -> f64 {
        let x = self.norm2() as f64;
        x.sqrt()
    }

    pub fn inf_norm(&self) -> i32 {
        let mut max = self.data[0];

        for i in 1..self.len() {
            if self.data[i] > max {
                max = self.data[i];
            }
        }

        max
    }
}

impl Matrix {
    pub fn from_slice(rows:usize, columns:usize, data:&[i32]) -> Matrix {
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
