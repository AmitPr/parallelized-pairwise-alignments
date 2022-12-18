use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    time::Instant,
};

use pairwise::{
    fast_lsa::FastLSAAligner, parallel::ParallelHirschbergAligner, sequential::HirschbergAligner,
    Aligner, Alignment,
};

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let file = args.get(1);
    let file = file.as_ref().map(|s| s.as_str()).unwrap_or_else(|| {
        eprintln!(
            "Usage: {} [file] [algorithm] [output_file (optional)]",
            args[0]
        );
        eprintln!("Available algorithms: all, hirschberg, parallel, fast_lsa");
        std::process::exit(1);
    });
    let (x, y) = read_sequences(file);

    let algorithm = args.get(2);
    let algorithm = algorithm.as_ref().map(|s| s.as_str()).unwrap_or_else(|| {
        eprintln!(
            "Usage: {} [file] [algorithm] [output_file (optional)]",
            args[0]
        );
        eprintln!("Available algorithms: all, hirschberg, parallel, fast_lsa");
        std::process::exit(1);
    });

    let output_file = args.get(3);
    let output_file = output_file
        .as_ref()
        .map(|s| std::fs::File::create(s).unwrap());

    match algorithm {
        "all" => {
            if output_file.is_some() {
                eprintln!("Output file not supported for running all algorithms.");
                std::process::exit(1);
            }
            run_all(&x, &y);
        }
        "hirschberg" => {
            run_hirschberg(&x, &y, output_file);
        }
        "parallel" => {
            run_parallel(&x, &y, output_file);
        }
        "fast_lsa" => {
            run_fast_lsa(&x, &y, output_file);
        }
        _ => {
            eprintln!("Unknown algorithm: {}", algorithm);
            eprintln!(
                "Usage: {} [file] [algorithm] [output_file (optional)]",
                args[0]
            );
            eprintln!("Available algorithms: all, hirschberg, parallel, fast_lsa");
            std::process::exit(1);
        }
    }
}

fn read_sequences(file: &str) -> (Vec<u8>, Vec<u8>) {
    // seq1 is first line, seq2 is second line
    let mut seq1 = String::new();
    let mut seq2 = String::new();

    let path = std::env::current_dir().unwrap().join(file);
    let file = File::open(path).unwrap_or_else(|_| panic!("Could not open file: {}", file));
    let mut reader = BufReader::new(file);
    reader
        .read_line(&mut seq1)
        .unwrap_or_else(|_| panic!("Could not read first line from file."));
    reader
        .read_line(&mut seq2)
        .unwrap_or_else(|_| panic!("Could not read second line from file."));
    (
        seq1.trim().as_bytes().to_vec(),
        seq2.trim().as_bytes().to_vec(),
    )
}

fn run_all(x: &[u8], y: &[u8]) {
    run_hirschberg(x, y, None);
    run_parallel(x, y, None);
    run_fast_lsa(x, y, None);
}

fn run_hirschberg(x: &[u8], y: &[u8], output_file: Option<std::fs::File>) {
    let score = |a: u8, b: u8| {
        if a == b'-' || b == b'-' {
            -2
        } else if a == b {
            2
        } else {
            -1
        }
    };
    let start = Instant::now();
    let aligner = HirschbergAligner::with(x, y, score);
    let aln = aligner.align();
    println!(
        "Hirschberg: Score: {} Time: {}",
        aln.score,
        start.elapsed().as_secs_f64(),
    );
    write_alignment(aln, output_file);
}

fn run_parallel(x: &[u8], y: &[u8], output_file: Option<std::fs::File>) {
    let score = |a: u8, b: u8| {
        if a == b'-' || b == b'-' {
            -2
        } else if a == b {
            2
        } else {
            -1
        }
    };
    let start = Instant::now();
    let aligner = ParallelHirschbergAligner::with(x, y, score);
    let aln = aligner.align();
    println!(
        "Parallel Hirschberg: Score: {} Time: {}",
        aln.score,
        start.elapsed().as_secs_f64(),
    );
    write_alignment(aln, output_file);
}

fn run_fast_lsa(x: &[u8], y: &[u8], output_file: Option<std::fs::File>) {
    let score = |a: u8, b: u8| {
        if a == b'-' || b == b'-' {
            -2
        } else if a == b {
            2
        } else {
            -1
        }
    };
    let start = Instant::now();
    let aligner = FastLSAAligner::with(x, y, score).with_block_size(500);
    let aln = aligner.align();
    println!(
        "Fast LSA: Score: {} Time: {}",
        aln.score,
        start.elapsed().as_secs_f64(),
    );
    write_alignment(aln, output_file);
}

fn write_alignment(aln: Alignment, output_file: Option<std::fs::File>) {
    if let Some(output_file) = output_file {
        let mut writer = std::io::BufWriter::new(output_file);
        let a = aln.alignment.iter().map(|(c1, c2)| *c1).collect::<String>() + "\n";
        let b = aln.alignment.iter().map(|(c1, c2)| c2).collect::<String>();
        writer
            .write_all(a.as_bytes())
            .unwrap_or_else(|_| panic!("Could not write to output file."));
        writer
            .write_all(b.as_bytes())
            .unwrap_or_else(|_| panic!("Could not write to output file."));
        println!("Wrote alignment to file.");
    }
}
