use std::fmt::Display;
use std::fmt;
use std::slice::Iter;

pub mod parser;
pub mod steinitz;
pub mod discrepancy;
mod graph;

pub type IntData = i32;
pub type Cost = i32;

#[derive(Hash, PartialEq, Eq, Clone)]
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
    delta_A: IntData,
    delta_b: IntData
}

pub enum ILPError {
    NoSolution,
    Unbounded,
    UnsupportedMatrix
}

impl ILP {
    pub fn new(mat:Matrix, b:Vector, c:Vector) -> Self {
        assert!(b.len() == mat.size.0);
        assert!(c.len() == mat.size.1);
        assert!(mat.size.0 > 0 && mat.size.1 > 0);

        let da = mat.max_abs_entry();
        let db = b.inf_norm();
    
        ILP {
            A: mat,
            b: b,
            c: c,
            delta_A: da,
            delta_b: db
        }
    }

    pub fn print_details(&self, prefix:&str) {
        println!("{}ILP details:", prefix);
        println!("{} -> constraints: {}",  prefix, self.A.size.0);
        println!("{} -> variables: {:3}",  prefix, self.A.size.1);
        println!("{} -> \u{0394}    = {}", prefix, self.delta_A);
        println!("{} -> \u{2016}b\u{2016}\u{221E} = {}", prefix, self.delta_b);
        if self.A.size.0 > 1 {
            println!("{} -> Matrix A:\n{}", prefix, self.A);
        } else {
            print!(  "{} -> Matrix A: {}",  prefix, self.A);
        }
        println!("{} -> b = {:?}",   prefix, self.b);
        println!("{} -> c = {:?}\n", prefix, self.c);
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

    pub fn unit(size:usize, dim:usize) -> Self {
        let mut data = vec![0; size];
        data[dim] = 1 as IntData;
        Vector {
            data: data
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

    pub fn max_distance(&self, v:&Vector, bound:IntData) -> bool {
        assert!(self.len() == v.len());

        for (&a,&b) in self.iter().zip(v.iter()) {
            if IntData::abs(a-b) > bound {
                return false;
            }
        }

        true
    }
}

impl fmt::Debug for Vector {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "{:?}", self.data)
    }
}

impl Matrix {
    pub fn zero(m:usize, n:usize) -> Self {
        let mut columns = Vec::with_capacity(n);

        for _ in 0..n {
            columns.push(Vector::zero(m));
        }

        Matrix {
            columns: columns,
            size: (m, n)
        }
    }

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

    pub fn num_cols(&self) -> usize {
        self.columns.len()
    }

    pub fn iter(&self) -> Iter<Vector> {
        self.columns.iter()
    }

    pub fn max_abs_entry(&self) -> IntData {
        self.iter().map(|col| col.inf_norm()).max().unwrap()
    }

    pub fn has_duplicate_columns(&self) -> bool {
        for (i,v) in self.iter().enumerate() {
            for c in self.iter().skip(i+1) {
                if v==c {
                    return true;
                }
            }
        }

        false
    }

    pub fn has_zero_columns(&self) -> bool {
        'column: for v in self.iter() {
            for &x in v.iter() {
                if x!=0 {
                    continue 'column;
                }
            }

            return true;
        }

        false
    }

    pub fn herdisc_upper_bound(&self) -> f32 {
        let (m,_) = self.size;
        let t = self.iter().map(|col| col.one_norm()).max().unwrap();

        let h = if m <= 699452 {
            2.0*f64::ln(2.0*m as f64)
        } else {
            5.32
        } as f32;

        let delta = self.max_abs_entry() as f32;

        f32::min(
            0.5 * h * f32::sqrt(m as f32) * delta,
            t as f32 // THM 7
        )
    }

    pub fn add_to_entry(&mut self, i:usize, j:usize, val:IntData) {
        self.columns[j].data[i] += val;
    }
}

impl Display for Matrix { 
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut str = "".to_string();
        let (m,n) = self.size;

        for i in 0..m {
            str.push_str("|");
            for j in 0..n {
                str.push_str(&format!(" {:3} ", self.columns[j].data[i]));
            }
            str.push_str("|\n");
        }

        write!(f, "{}", str)
    }
}
