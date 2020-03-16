pub mod ilp;
use ilp::*;
use clap::{App, Arg};

fn main() {
    let matches = App::new("IntOpt ILP Solver")
        .version(env!("CARGO_PKG_VERSION"))
        .author("Tim Weißenfels <tim.weissenfels@stu.uni-kiel.de>")
        .about("Solves integer linear programs (ILPs).")
        .arg(
            Arg::with_name("algorithm")
                .short("a")
                .long("algorithm")
                .value_name("ALGORITHM")
                .default_value("steinitz")
                .possible_values(&["steinitz", "discrepancy"])
                .help("Sets the algorithm to solve the ILP.")
                .takes_value(true),
        )
        .get_matches();

    /*let mat = Matrix::from_slice(6, 6, &[2,0,0,0,0,0, 0,3,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,1,1,0, 1,0,0,0,0,1]);
    let b = Vector::from_slice(&[3, 0, 0, 1, 1, 1]);
    let c = Vector::from_slice(&[1, 2, 3, 4, 3, 2]);*/
    let mat = Matrix::from_slice(3, 3, &[1,0,0,0,2,0,0,0,1]);
    let b = Vector::from_slice(&[5, 6, 5]);
    let c = Vector::from_slice(&[1, 2, 3]);

    let ilp = ILP::new(mat, b, c);
    let res = match matches.value_of("algorithm") {
        Some("steinitz") => steinitz::solve(&ilp),
        Some("discrepancy") => discrepancy::solve(&ilp),
        _ => panic!()
    };

    match res {
        Ok(x) => println!("Found a solution! x={:?}", x),
        Err(ILPError::NoSolution) => println!("The ILP has no solution."),
        Err(ILPError::Unbounded)  => println!("The ILP is unbounded."),
        Err(_) => println!("This ILP could not be solved.")
    }
}
