# Parallelized Pairwise Alignments
This project contains the implementations of three pairwise alignment algorithms:
* Hirschberg's algorithm
* Parallelized Hirschberg's algorithm
* FastLSA algorithm

## Requirements
* Rust (https://www.rust-lang.org/tools/install)
Please install the latest stable version of Rust.

## Installation
First clone this repository:
```
git clone https://github.com/AmitPr/parallelized-pairwise-alignments.git
```
Then, install the project executable:
```
cd parallelized-pairwise-alignments
cargo install --path .
```
This will install the executable `pairwise` in your `~/.cargo/bin` directory, which should be in your `PATH` by default.

## Usage
```
pairwise [input_file] [algorithm] [output_file (optional)]
```
* `input_file`: The path to the input file. The input file should contain two sequences, each on a separate line. These are the sequences that will be aligned.
* `algorithm`: The name of the algorithm to use. The available algorithms are: `hirschberg`, `parallel`, `fast_lsa`, and `all` for running all algorithms.
* `output_file`: The path to the output file. If not specified, there will be no output and only statistics will be printed.

If using the `all` algorithm, please do not specify an output file, as only the statistics will be printed. The `all` option can be used to compare the performance of the different algorithms.

## Example
```
cd parallelized-pairwise-alignments
pairwise benches/short.txt fast_lsa output.txt
cat output.txt
# Output will be:
AGTACGCA
--TATGC-
```
