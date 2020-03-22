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
enum Constraint {
    Equation   { left: Sum, right: Sum },
    Inequality { left: Sum, right: Sum, leq:bool }
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
    let objective_tree;
    let constraints_tree;

    {
        let mut iterator = file.into_inner();
        maximize = iterator.next().unwrap().as_str().to_lowercase() == "maximize";
        objective_tree = iterator.next().unwrap();
        constraints_tree = iterator.next().unwrap();
    }

    // find variables
    {
        let vars1 = find_variables(&objective_tree);
        let vars2 = find_variables(&constraints_tree);
        let mut list = Vec::new();
        
        for var in vars1.iter().chain(vars2.iter()) {
            if !variables.contains_key(var) {
                list.push(var);
                variables.insert(var.clone(), variables.len());
            }
        }

        println!("Variables: {:?}", list);
    }

    let constraints = get_constraints(constraints_tree);
    let inequalities = constraints.iter().filter(|c| matches!(c, Constraint::Inequality{..})).count();
    let m = constraints.len();
    let n = variables.len() + inequalities; // a slack var for every inequality
    let mut a = Matrix::zero(m, n);
    let mut b = Vector::zero(m);
    let mut c = Vector::zero(n);

    // objective -> c Vector
    for m in multiple_sum(objective_tree).1 {
        let i = *variables.get(&m.1).unwrap();
        if maximize {
            c.data[i] += m.0;
        } else {
            c.data[i] -= m.0;
        }
        
    }

    // constraints -> A matrix
    let mut slack = 0;
    for (row, c) in constraints.iter().enumerate() {
        let (left, right) = match c {
            Constraint::Equation{ left, right } => (left, right),
            Constraint::Inequality{ left, right, leq } => {
                let j = variables.len() + slack;
                slack += 1;
                a.add_to_entry(row, j, if *leq {1} else {-1});
                (left,  right)
            }
        };

        b.data[row] = right.0 - left.0;
        for m in left.1.iter() {
            let j = *variables.get(&m.1).unwrap();
            a.add_to_entry(row, j, m.0);
        }
        for m in right.1.iter() {
            let j = *variables.get(&m.1).unwrap();
            a.add_to_entry(row, j, -m.0);
        }
    }

    if slack > 0 {
        println!("Introduced {} slack variables.", slack);
    }

    println!();

    Ok(ILP::with_named_vars(a,b,c,variables.drain().collect()))
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

fn constraint(pair: Pair<Rule>) -> Constraint {
    let rule = pair.as_rule();
    let mut iter = pair.into_inner();
    let left  = multiple_sum(iter.next().unwrap());
    let right = multiple_sum(iter.next().unwrap());

    match rule {
        Rule::equation => Constraint::Equation { left: left, right: right },
        Rule::leq      => Constraint::Inequality { left: left, right: right, leq: true },
        Rule::geq      => Constraint::Inequality { left: left, right: right, leq: false },
        _              => unreachable!()
    }
}

fn get_constraints(pair: Pair<Rule>) -> Vec<Constraint> {
    assert_eq!(pair.as_rule(), Rule::constraints);

    fn f(v:&mut Vec<Constraint>, pair:Pair<Rule>) {
        for p in pair.into_inner() {
            match p.as_rule() {
                Rule::equation    => v.push(constraint(p)),
                Rule::leq         => v.push(constraint(p)),
                Rule::geq         => v.push(constraint(p)),
                Rule::constraints => f(v, p),
                _                 => unreachable!()
            }
        }
    }

    let mut v = Vec::new();
    f(&mut v, pair);
    v
}
