static SCORER: fn(u8, u8) -> i32 = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pairwise::fast_lsa::FastLSAAligner;
use pairwise::parallel::ParallelHirschbergAligner;
use pairwise::sequential::HirschbergAligner;
use pairwise::Aligner;

fn bench_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("Alignment");
    for (file, bs) in &[("short.txt", 10), ("medium.txt", 50), ("long.txt", 500)] {
        let (seq_1, seq_2) = read_sequences(file);
        group.bench_with_input(
            format!("sequential/{file}"),
            &(seq_1.clone(), seq_2.clone()),
            |b, (seq_1, seq_2)| {
                b.iter(|| {
                    let aligner = HirschbergAligner::with(
                        black_box(seq_1),
                        black_box(seq_2),
                        black_box(SCORER),
                    );
                    aligner.align();
                })
            },
        );
        group.bench_with_input(
            format!("parallel/{file}"),
            &(seq_1.clone(), seq_2.clone()),
            |b, (seq_1, seq_2)| {
                b.iter(|| {
                    let aligner = ParallelHirschbergAligner::with(
                        black_box(seq_1),
                        black_box(seq_2),
                        black_box(SCORER),
                    );
                    aligner.align();
                })
            },
        );
        group.bench_with_input(
            format!("fastlsa/{file}"),
            &(seq_1.clone(), seq_2.clone()),
            |b, (seq_1, seq_2)| {
                b.iter(|| {
                    let aligner =
                        FastLSAAligner::with(black_box(seq_1), black_box(seq_2), black_box(SCORER))
                            .with_block_size(*bs);
                    aligner.align();
                })
            },
        );
    }
}

fn read_sequences(file: &str) -> (Vec<u8>, Vec<u8>) {
    // seq1 is first line, seq2 is second line
    let mut seq1 = String::new();
    let mut seq2 = String::new();
    // cwd + "/benches/{file}"
    let path = std::env::current_dir().unwrap().join("benches").join(file);
    let file = File::open(path).unwrap();
    let mut reader = BufReader::new(file);
    reader.read_line(&mut seq1).unwrap();
    reader.read_line(&mut seq2).unwrap();
    (
        seq1.trim().as_bytes().to_vec(),
        seq2.trim().as_bytes().to_vec(),
    )
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_group
}

criterion_main!(benches);
