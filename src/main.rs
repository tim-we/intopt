extern crate pest;
#[macro_use] extern crate pest_derive;
#[macro_use] extern crate matches;

pub mod ilp;
use ilp::*;
use clap::{App, Arg};

fn main() {
    let matches = App::new("IntOpt ILP Solver")
        .version(env!("CARGO_PKG_VERSION"))
        .author("Tim Weißenfels <tim.weissenfels@stu.uni-kiel.de>")
        //.about("Solves integer linear programs (ILPs).\nhttps://github.com/tim-we/intopt")
        .about(
            &format!("{}\n{}\n{}",
                env!("CARGO_PKG_DESCRIPTION"),
                "max { <c,x> | Ax=b, 0\u{2264}x, x\u{2208}\u{2124}\u{207F} }",
                env!("CARGO_PKG_REPOSITORY")
            )[..])
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
        .arg(
            Arg::with_name("input")
                .takes_value(true)
                .value_name("FILE")
                .help(".ilp input file")
                .required(true)
        )
        .get_matches();

    let ilp = parser::parse_file(matches.value_of("input").unwrap()).unwrap();
    ilp.print_details("");

    let res = match matches.value_of("algorithm") {
        Some("steinitz")    => steinitz::solve(&ilp),
        Some("discrepancy") => discrepancy::solve(&ilp),
        _ => panic!()
    };

    println!();

    match res {
        Ok(x) => {
            println!("Solution:");
            ilp.print_solution(&x)
        },
        Err(ILPError::NoSolution) => println!("The ILP has no solution."),
        Err(ILPError::Unbounded)  => println!("The ILP is unbounded."),
        Err(_) => println!("This ILP could not be solved.")
    }
}
