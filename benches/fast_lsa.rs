static SCORER: fn(u8, u8) -> i32 = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pairwise::fast_lsa::FastLSAAligner;
use pairwise::Aligner;
use pprof::criterion::{Output, PProfProfiler};

fn bench_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("FastLSA");
    for file in &["short.txt", "medium.txt"] {
        for block_size in &[50, 75, 150, 250, 400] {
            let (seq_1, seq_2) = read_sequences(file);
            group.bench_with_input(
                format!("{block_size}-{file}"),
                &(seq_1, seq_2),
                |b, (seq_1, seq_2)| {
                    b.iter(|| {
                        let aligner = FastLSAAligner::with(
                            black_box(seq_1),
                            black_box(seq_2),
                            black_box(SCORER),
                        )
                        .with_block_size(*block_size);
                        aligner.align();
                    })
                },
            );
        }
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
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    targets = bench_group
}

criterion_main!(benches);
