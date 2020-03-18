use pest::Parser;
use pest::iterators::Pair;
use std::fs;
use super::{ILP, Vector, Matrix};

type Map<K,V> = hashbrown::HashMap<K,V>;
type Set<T> = hashbrown::HashSet<T>;

#[derive(Parser)]
#[grammar = "ilp.pest"]
pub struct ILPFileParser;

struct Multiple(i32,String);
struct Sum(i32,Vec<Multiple>);
struct Equation {
    left:  Sum,
    right: Sum
}

pub fn parse_file(file:&str) -> Result<ILP, ()> {
    println!("Reading file {}...", file);
    let unparsed_file = fs::read_to_string(file).expect("cannot read file");

    // parse file
    println!("Parsing file...");
    let file = ILPFileParser::parse(Rule::ilp, &unparsed_file)
                .expect("unsuccessful parse")
                .next().unwrap();
    
    let mut variables = Map::<String, usize>::new();
    let maximize;
    let objective;
    let constraints;

    {
        let mut iterator = file.into_inner();
        maximize = iterator.next().unwrap().as_str().to_lowercase() == "maximize";
        objective = iterator.next().unwrap();
        constraints = iterator.next().unwrap();
    }

    // find variables
    {
        let vars1 = find_variables(&objective);
        let vars2 = find_variables(&constraints);
        let mut list = Vec::new();
        
        for var in vars1.iter().chain(vars2.iter()) {
            if !variables.contains_key(var) {
                list.push(var);
                variables.insert(var.clone(), variables.len());
            }
        }

        println!("Variables: {:?}", list);
    }

    let equations = get_equations(constraints);
    let m = equations.len();
    let n = variables.len();
    let mut a = Matrix::zero(m, n);
    let mut b = Vector::zero(m);
    let mut c = Vector::zero(n);

    // objective -> c Vector
    for m in multiple_sum(objective).1 {
        let i = *variables.get(&m.1).unwrap();
        if maximize {
            c.data[i] += m.0;
        } else {
            c.data[i] -= m.0;
        }
        
    }

    // constraints -> A matrix
    for (row, eq) in equations.iter().enumerate() {
        b.data[row] = eq.right.0 - eq.left.0;
        for m in eq.left.1.iter() {
            let j = *variables.get(&m.1).unwrap();
            a.add_to_entry(row, j, m.0);
        }
        for m in eq.right.1.iter() {
            let j = *variables.get(&m.1).unwrap();
            a.add_to_entry(row, j, -m.0);
        }
    }

    println!("Parsing successful.\n");

    Ok(ILP::new(a,b,c))
}

fn find_variables(tree: &Pair<Rule>) -> Vec<String> {
    let mut set = Set::<String>::new();
    let mut list = Vec::new();

    let var_strs = tree.clone()
        .into_inner()
        .flatten()
        .filter(|p| p.as_rule() == Rule::variable)
        .map(|p| p.as_str());

    for str in var_strs {
        let s = String::from(str);
        if !set.contains(&s) {
            set.insert(s.clone());
            list.push(s);
        }
    }

    list
}

fn multiple_sum(pair: Pair<Rule>) -> Sum {
    assert_eq!(pair.as_rule(), Rule::sum);

    fn build_sum(sum:&mut Sum, pair: Pair<Rule>) {
        for p in pair.into_inner() {
            match p.as_rule() {
                Rule::integer  => sum.0 += p.as_str().parse::<i32>().unwrap(),
                Rule::multiple => sum.1.push(multiple(p)),
                Rule::term     => build_sum(sum, p),
                Rule::sum      => build_sum(sum, p),
                _              => unreachable!()
            }
        }
    }

    let mut sum = Sum(0, Vec::new());
    build_sum(&mut sum, pair);
    sum
}

fn multiple(pair: Pair<Rule>) -> Multiple {
    assert_eq!(pair.as_rule(), Rule::multiple);

    let mut var_name = "".to_string();
    let mut multiple = 1;

    for p in pair.into_inner() {
        match p.as_rule() {
            Rule::integer  => multiple = p.as_str().parse().unwrap(),
            Rule::variable => var_name = p.as_str().to_string(),
            _ => unreachable!()
        }
    }

    Multiple(multiple, var_name)
}

fn equation(pair: Pair<Rule>) -> Equation {
    assert_eq!(pair.as_rule(), Rule::equation);

    let mut iter = pair.into_inner();
    let left = iter.next().unwrap();
    let right = iter.next().unwrap();

    Equation {
        left: multiple_sum(left),
        right: multiple_sum(right)
    }
}

fn get_equations(pair: Pair<Rule>) -> Vec<Equation> {
    assert_eq!(pair.as_rule(), Rule::constraints);

    fn f(v:&mut Vec<Equation>, pair:Pair<Rule>) {
        for p in pair.into_inner() {
            match p.as_rule() {
                Rule::equation    => v.push(equation(p)),
                Rule::constraints => f(v, p),
                _                 => unreachable!()
            }
        }
    }

    let mut v = Vec::new();
    f(&mut v, pair);
    v
}
