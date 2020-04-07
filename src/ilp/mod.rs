use std::fmt::Display;
use std::fmt;
use std::slice::Iter;

pub mod parser;
pub mod steinitz;
pub mod discrepancy;
mod graph;

pub type IntData = i32;
pub type Cost = i32;
pub type VarMapping = (String, usize);

#[derive(Hash, PartialEq, Eq, Clone)]
pub struct Vector {
    data: Vec<IntData>
}

#[derive(Clone)]
pub struct Matrix {
    columns: Vec<Vector>,
    size: (usize, usize) // rows, columns or (m,n)
}

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct ILP {
    pub A: Matrix,
    pub b: Vector,
    pub c: Vector,
    delta_A: IntData,
    delta_b: IntData,
    named_variables: Vec<VarMapping>
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

        let da = mat.max_abs_entry();
        let db = b.inf_norm();

        assert!(da >= 0 && db >= 0);
    
        ILP {
            A: mat,
            b: b,
            c: c,
            delta_A: da,
            delta_b: db,
            named_variables: Vec::new()
        }
    }

    pub fn with_named_vars(mat:Matrix, b:Vector, c:Vector, vars:Vec<VarMapping>) -> Self {
        let mut ilp = ILP::new(mat, b, c);

        for (s, idx) in vars.iter() {
            assert!(s.len() > 0);
            assert!(idx < &ilp.c.len());
        }

        let mut variables = vars;
        variables.sort_by(|a,b| a.1.cmp(&b.1));

        ilp.named_variables = variables;
        ilp
    }

    pub fn print_details(&self) {
        println!("ILP details:");
        println!(" -> constraints: {}", self.A.size.0);
        println!(" -> variables: {:3}", self.A.size.1);
        let list:Vec<&String> = self.named_variables.iter().map(|(s,_)| s).collect();
        print!("    {:?}", list);
        let slacks = self.A.size.1 - list.len();
        if slacks > 0 {
            println!(" + {} slack variables", slacks);
        } else {
            println!();
        }
        println!(" -> \u{0394}    = {}", self.delta_A);
        println!(" -> \u{2016}b\u{2016}\u{221E} = {}", self.delta_b);
        if self.A.size.0 > 1 {
            println!(" -> Matrix A:\n{}", self.A);
        } else {
            print!(  " -> Matrix A: {}", self.A);
        }
        println!(" -> b = {:?}", self.b);
        println!(" -> c = {:?}\n", self.c);
    }

    pub fn print_solution(&self, x:&Vector) {
        if self.named_variables.len() == 0 {
            println!(" x={:?}", x);
        } else {
            for (name, idx) in self.named_variables.iter() {
                println!(" {} = {}", name, x.data[*idx]);
            }
        }
    }

    pub fn simplify(self) -> Self {
        assert!(self.A.columns.len() > 1);
        
        let mut mat = Matrix {
            columns: Vec::with_capacity(self.A.size.1 - 1),
            size: (self.b.len(), 0)
        };
    
        let mut c = Vector {
            data: Vec::new()
        };
        
        let mut var_names:Vec<Option<String>> = vec![None; self.A.size.1];
        self.named_variables.iter().for_each(|(str, i)| var_names[*i] = Some(str.clone()));
        
        let mut skip = Vec::new();
        for (i, col1) in self.A.iter().enumerate() {
            if skip.contains(&i) {
                continue;
            }
    
            let mut best = (col1, self.c.data[i]);
            for (j, col2) in self.A.iter().enumerate().skip(i+1) {
                if col1 == col2 {
                    let cost = self.c.data[j];
                    
                    // keep column with highest cost/weight
                    let removed = if cost > best.1 {
                        best = (col2, cost);
                        var_names.remove(i)
                    } else {
                        var_names.remove(j)
                    };

                    if let Some(name) = removed {
                        println!("    {} = 0", name);
                    }

                    skip.push(j);
                }
            }
            
            mat.columns.push(best.0.clone());
            mat.size.1 += 1;
            c.data.push(best.1);
        }

        let mappings = var_names.into_iter()
            .enumerate()
            .map(|(i, o)| match o {
                Some(str) => Some((str, i)),
                None      => None
            })
            .flatten()
            .collect();

        println!(" -> Removed {} column(s).", skip.len());
    
        ILP::with_named_vars(mat, self.b.clone(), c, mappings)
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
        debug_assert!(self.len() == other.len());
        let mut v = Vec::with_capacity(self.len());

        for (x1,x2) in self.iter().zip(other.iter()) {
            v.push(x1 + x2);
        }

        Vector {
            data: v
        }
    }

    pub fn dot(&self, other: &Vector) -> IntData {
        debug_assert!(self.len() == other.len());
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
        debug_assert!(self.len() == v.len());

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

    pub fn non_negative(&self) -> bool {
        for c in self.columns.iter() {
            if c.iter().filter(|&&x| x < 0).count() > 0 {
                return false;
            }
        }

        true
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
