use std::sync::{Arc, Mutex};

use crate::{Aligner, Alignment};

const DEFAULT_BLOCKSIZE: usize = 150;

type Cache<T> = Arc<Vec<Mutex<Vec<T>>>>;

pub struct FastLSAAligner<'a> {
    a: &'a [u8],
    b: &'a [u8],
    score: fn(u8, u8) -> i32,
    block_size: usize,
}

impl<'a> Aligner<'a> for FastLSAAligner<'a> {
    fn with(a: &'a [u8], b: &'a [u8], score: fn(u8, u8) -> i32) -> Self {
        FastLSAAligner {
            a,
            b,
            score,
            block_size: DEFAULT_BLOCKSIZE,
        }
    }

    fn align(&self) -> Alignment {
        let score = self.score;
        // check if alignment fits within block size
        if self.a.len() <= self.block_size && self.b.len() <= self.block_size {
            return crate::sequential::HirschbergAligner::with(self.a, self.b, score).align();
        }

        let num_grid_rows = (self.b.len() + self.block_size - 1) / self.block_size;
        let num_grid_cols = (self.a.len() + self.block_size - 1) / self.block_size;
        let num_cached_rows = num_grid_rows + 1;
        let num_cached_cols = num_grid_cols + 1;

        let row_cache: Cache<i32> = Arc::new(
            std::iter::repeat_with(|| Mutex::new(vec![0; self.a.len() + 1]))
                .take(num_cached_rows)
                .collect::<Vec<_>>(),
        );
        let col_cache: Cache<i32> = Arc::new(
            std::iter::repeat_with(|| Mutex::new(vec![0; self.b.len() + 1]))
                .take(num_cached_cols)
                .collect::<Vec<_>>(),
        );

        let visited_blocks = Arc::new(Mutex::new(vec![
            vec![false; num_cached_cols];
            num_cached_rows
        ]));

        let done_blocks = Arc::new(Mutex::new(vec![
            vec![false; num_cached_cols];
            num_cached_rows
        ]));

        /*
         * 1. Fill in the cache grid
         */

        // First row and first column are always just insertion and deletion scores
        {
            let mut row_cache = row_cache[0].lock().unwrap();
            for i in 1..=self.a.len() {
                row_cache[i] = row_cache[i - 1] + score(self.a[i - 1], b'-');
            }
        }
        {
            let mut col_cache = col_cache[0].lock().unwrap();
            for i in 1..=self.b.len() {
                col_cache[i] = col_cache[i - 1] + score(b'-', self.b[i - 1]);
            }
        }

        // Fill in the rest of the cache
        self.forward(
            (0, 0),
            row_cache.clone(),
            col_cache.clone(),
            visited_blocks,
            done_blocks,
        );

        /*
         * 2. Backtrace
         */
        let mut last_a = &self.a[(num_grid_cols - 1) * self.block_size..];
        let mut last_b = &self.b[(num_grid_rows - 1) * self.block_size..];
        let mut start_row =
            row_cache[num_grid_rows - 1].lock().unwrap()[self.a.len() - last_a.len()..].to_vec();
        let mut start_col =
            col_cache[num_grid_cols - 1].lock().unwrap()[self.b.len() - last_b.len()..].to_vec();
        let mut start_idx = (last_a.len(), last_b.len());
        let mut next_row = num_grid_rows - 1;
        let mut next_col = num_grid_cols - 1;

        let mut alignment = vec![];
        let alignment_score = *row_cache.last().unwrap().lock().unwrap().last().unwrap();
        loop {
            let result = Self::backtrace_grid_cell(
                (next_row, next_col),
                last_a,
                last_b,
                start_row,
                start_col,
                start_idx,
                score,
            );
            alignment.splice(0..0, result.0.alignment);

            if next_row == 0 && next_col == 0 {
                break;
            }
            // (sub_alignment, start_idx, next_row, next_col)
            start_idx = result.1;
            if next_col != 0 && start_idx.0 == 0{
                start_idx.0 = self.block_size
            }

            if next_row != 0 && start_idx.1 == 0{
                start_idx.1 = self.block_size
            }

            next_row = result.2;
            next_col = result.3;

            last_a = &self.a
                [next_col * self.block_size..((next_col + 1) * self.block_size).min(self.a.len())];
            last_b = &self.b
                [next_row * self.block_size..((next_row + 1) * self.block_size).min(self.b.len())];

            start_row = row_cache[next_row].lock().unwrap()[next_col * self.block_size
                ..((next_col + 1) * self.block_size + 1).min(self.a.len() + 1)]
                .to_vec();

            start_col = col_cache[next_col].lock().unwrap()[next_row * self.block_size
                ..((next_row + 1) * self.block_size + 1).min(self.b.len() + 1)]
                .to_vec();
        }

        Alignment {
            alignment,
            score: alignment_score,
        }
    }
}

impl<'a> FastLSAAligner<'a> {
    pub fn with_block_size(mut self, block_size: usize) -> Self {
        self.block_size = block_size;
        self
    }

    pub fn forward(
        &self,
        cell: (usize, usize),
        row_cache: Cache<i32>,
        col_cache: Cache<i32>,
        visited_blocks: Arc<Mutex<Vec<Vec<bool>>>>,
        done_blocks: Arc<Mutex<Vec<Vec<bool>>>>,
    ) {
        let (row, col) = cell;
        let score = self.score;
        let num_grid_rows = (self.b.len() + self.block_size - 1) / self.block_size;
        let num_grid_cols = (self.a.len() + self.block_size - 1) / self.block_size;
        assert!(row < num_grid_rows, "row out of bounds");
        assert!(col < num_grid_cols, "col out of bounds");
        visited_blocks.lock().unwrap()[row][col] = true;

        let row_start = col * self.block_size;
        let row_end = ((col + 1) * self.block_size).min(self.a.len());
        let col_start = row * self.block_size;
        let col_end = ((row + 1) * self.block_size).min(self.b.len());

        // println!("Filling in cell ({}, {})", row, col);
        // println!(
        //     "rs: {}, re: {}, cs: {}, ce: {}",
        //     row_start, row_end, col_start, col_end
        // );

        let start_row = row_cache[row].lock().unwrap()[row_start..=row_end].to_vec();
        let start_col = col_cache[col].lock().unwrap()[col_start..=col_end].to_vec();
        let mut output_row = vec![0; row_end - row_start + 1];
        let mut output_col = vec![0; col_end - col_start + 1];
        Self::fill_grid_block(
            &self.a[row_start..row_end],
            &self.b[col_start..col_end],
            start_row,
            start_col,
            output_row.as_mut_slice(),
            output_col.as_mut_slice(),
            score,
        );

        row_cache[row + 1].lock().unwrap()[row_start..=row_end].copy_from_slice(&output_row);
        col_cache[col + 1].lock().unwrap()[col_start..=col_end].copy_from_slice(&output_col);

        let south_ready: bool;
        let east_ready: bool;
        {
            let mut done_blocks = done_blocks.lock().unwrap();
            let visited_blocks = visited_blocks.lock().unwrap();
            south_ready = row + 1 < num_grid_rows
                && !visited_blocks[row + 1][col]
                && (col == 0 || done_blocks[row + 1][col - 1]);
            east_ready = col + 1 < num_grid_cols
                && !visited_blocks[row][col + 1]
                && (row == 0 || done_blocks[row - 1][col + 1]);

            done_blocks[row][col] = true;
        }

        rayon::join(
            || {
                if south_ready {
                    self.forward(
                        (row + 1, col),
                        row_cache.clone(),
                        col_cache.clone(),
                        visited_blocks.clone(),
                        done_blocks.clone(),
                    );
                }
            },
            || {
                if east_ready {
                    self.forward(
                        (row, col + 1),
                        row_cache.clone(),
                        col_cache.clone(),
                        visited_blocks.clone(),
                        done_blocks.clone(),
                    );
                }
            },
        );
    }

    pub fn fill_grid_block(
        a: &[u8],
        b: &[u8],
        start_row: impl AsRef<[i32]>,
        start_col: impl AsRef<[i32]>,
        mut output_row: impl AsMut<[i32]>,
        mut output_col: impl AsMut<[i32]>,
        score: fn(u8, u8) -> i32,
    ) {
        let n = a.len();
        let m = b.len();
        // first index already computed in start_{row,col}
        let mut scores = (vec![0; n + 1], vec![0; n + 1]);
        scores.0 = start_row.as_ref().to_vec();
        output_col.as_mut()[0] = start_row.as_ref()[n];
        for i in 0..m {
            scores.1[0] = start_col.as_ref()[i + 1];

            for j in 0..n {
                scores.1[j + 1] = {
                    let subst = scores.0[j] + score(a[j], b[i]);
                    let del = scores.1[j] + score(a[j], b'-');
                    let ins = scores.0[j + 1] + score(b'-', b[i]);
                    subst.max(del).max(ins)
                };
                // Fill output column
                if j == n - 1 {
                    output_col.as_mut()[i + 1] = scores.1[j + 1];
                }
            }
            scores.0 = scores.1.clone();
        }
        // Fill output row
        output_row.as_mut().copy_from_slice(&scores.1);
    }

    /// Returns: (alignment, next_start_idx, next_row, next_col)
    pub fn backtrace_grid_cell(
        cell: (usize, usize),
        a: &[u8],
        b: &[u8],
        start_row: impl AsRef<[i32]>,
        start_col: impl AsRef<[i32]>,
        start_idx: (usize, usize),
        score: fn(u8, u8) -> i32,
    ) -> (Alignment, (usize, usize), usize, usize) {
        #[derive(Debug, Clone, Copy, PartialEq)]
        enum Direction {
            None,
            Diagonal,
            Left,
            Up,
        }
        let (row, col) = cell;
        // println!(
        //     "Backward: row: {}, col: {}, start: ({}, {}), a: {:?}, b: {:?}",
        //     row,
        //     col,
        //     start_idx.0,
        //     start_idx.1,
        //     String::from_utf8_lossy(a),
        //     String::from_utf8_lossy(b)
        // );
        let n = a.len();
        let m = b.len();
        let mut scores = vec![vec![0; n + 1]; m + 1];
        let mut backtrack = vec![vec![Direction::None; n + 1]; m + 1];
        if row == 0 {
            for j in 0..n {
                backtrack[0][j + 1] = Direction::Left;
            }
        }
        if col == 0 {
            for i in 0..m {
                backtrack[i + 1][0] = Direction::Up;
            }
        }
        scores[0] = start_row.as_ref().to_vec();
        for i in 0..m {
            scores[i + 1][0] = start_col.as_ref()[i + 1];
            for j in 0..n {
                let nexts = vec![
                    (scores[i][j] + score(a[j], b[i]), Direction::Diagonal), // substitution
                    (scores[i][j + 1] + score(b'-', b[i]), Direction::Up),   // deletion
                    (scores[i + 1][j] + score(a[j], b'-'), Direction::Left), // insertion
                ];
                let max_score = *nexts.iter().map(|(score, _)| score).max().unwrap();
                backtrack[i + 1][j + 1] = nexts
                    .iter()
                    .find(|(score, _)| *score == max_score)
                    .unwrap()
                    .1;
                scores[i + 1][j + 1] = max_score;
            }
        }

        let mut i = start_idx.0;
        let mut j = start_idx.1;
        let mut alignment = Vec::new();
        let condition = |i: usize, j: usize| {
            if row == 0 && col == 0 {
                i > 0 || j > 0
            } else {
                i > 0 && j > 0
            }
        };
        while backtrack[j][i] != Direction::None {
            match backtrack[j][i] {
                Direction::Diagonal => {
                    alignment.push((a[i - 1] as char, b[j - 1] as char));
                    i -= 1;
                    j -= 1;
                }
                Direction::Left => {
                    alignment.push((a[i - 1] as char, '-'));
                    i -= 1;
                }
                Direction::Up => {
                    alignment.push(('-', b[j - 1] as char));
                    j -= 1;
                }
                _ => unreachable!(),
            }
        }
        alignment.reverse();
        let alignment = Alignment {
            score: scores[start_idx.1][start_idx.0],
            alignment,
        };
        if row == 0 && col == 0 {
            (alignment, (0, 0), 0, 0)
        } else {
            let next_row = if j > 0 || row == 0 { row } else { row - 1 };
            let next_col = if i > 0 || col == 0 { col } else { col - 1 };
            (alignment, (i, j), next_row, next_col)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Aligner;

    #[test]
    fn doit() {
        let a = b"ATCTAACTATTCCCTGTGCCTTATGGGGGCCTGCGCTATCTGCCTGTCGAACCATAGGACTCGCGCCAGCGCGCAGGCTTGGATCGAGGTGAAATCTCCGGGGCCTAAGACCACGAGCGTCTGGCGTCTTGGCTAACCCCCCTACATGCTGTTATAGACAATCAGTGGAAACCCGGTGCCAGGGGGTGGAGTGACCTTAAGTCAGGGACGATATTAATCGGAAGGAGTATTCAACGCAATGAAGCCGCAGGGTTGGCGTGGGAATGGTGCTTCTGTCCAAGCAGGTAAGGGCATGAGGCCGCAACCGTCCCCCAAGCGTACAGGGTGCACTTTGCAACGATTTCGGAGTCCGGAGACTCGCTGTTTTCGAAATTTGCGCTCAAGGGCGGGTATTGAACCAGGCTTACGCCCAAGAACGTAGCAAGGTGACTCAAACAAGGTACATCTTGCCCGCGTTTCACACGAATCAAGTTGGAGGTTATGGAGCATAGTAACACGTGGGCGGCCAGTGGTCGGTTGCTACACCCCTGCCGCAACGTTGAAGGTCCCGGATTAGACTGGCTGGACCCATGCCGTGACACCCGTCACACTCCATTACCGTCTGCGGGTCACGGCTTGTTGTGGACTGGATTGCCATTCTCTCAGTGTATTACGCAGGCCGGCGCGCGGGTCCCATGTAAACCTGTCATAGCTTACCTGACTCTACTTGGAAGTGTGGCTAGGCCTTTGCCCACGCACCTGGTCGGTCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCGCTGCAAGAATCTCTACCTGCTTTACAAAGCGCTGGGTCCTACTCCAGCGGGATGTTTTATCTAAACACGATGAGAGGAGTATTCGTCAGGCCACATGGCTTTCTTGTCCTGGTCGGATCCATCGTTGGCGCCCGACCCCCCCACTCCGTAGTGAGTTCTTCGTCCGAGCCATTGCATGCCAGATCGGCAGACAGATAGCGGATCCAGTATATCCCTGGAAGCTATAGACGCACAGGTTGGAATCCTAAGCGAAGTCGCGCGTCCGAACCCAGCTCTACTTTAGTGGCCACGGGTTCTGGTCCCCCCGGGCCGCGGAACCGATTAGGGCCATGTACAACAATACTTATTAGTCACCTTTCAGACACGATCTCCCTGCTCAGTGGTATATGGTTCCTGCTATAATTAGCCACCCTCATAAGTTGCACTACTTCTGCGACCCAAGTGCACCCTTACCACGAAGACAGGATTGTCCGATCCCATACTGCGGCCTTGGCAGGGGGTTCGCAAGTCCCACCCCAAACGATGCTGAAGGCTCAGGTTACACAGGCACAAGTGCTATATACGCGAGTTCCCGCTCTTAACCTGGACCGAATGCGGGATCATGCATCGTACCACTGTGTTCGTGTCATCTAGGACGGGCGCAAAGGATACATAGTTCAATCAAGAATACCTTGTATTATTGTACACCTACCGGTCACCAGCCAACAATGTGCGGACGGCGTTGCGACTTGCTGGGCCTGATCTCACCGCCCTAGATACCGCACACTGGGCAATACGAGGTAAAGCCAGTCACCCAGTGTCGATCAACAGCTGACGTAACGGTAAGAGGCTCACAAAATCGCACCGCCGGCGTCCCCTGGGTATTTTACGTCAGCATCGGGTGGACTGGCATGAATCTTTACTCCCAGGCGGAAACGGGTGCGTGGACAAGCGAGCAGCAAACGAAAATTCCTGGCCTGCTTGGTGTCTCGTATCCCTCTTGGAGATCGAGGAAATGTTTCACGACCAAGGGAAAGGTCGCCCTACGAAATAGATTTGCGCTACTGTCCGCATAAGGAGTCCGGTGTAGCGAAGGATGAAGGCGACCCTAGGTAGCAACCGCCGGCTTCGGCGGTAAGGTATCACTCAGGAAGCAGGCACGGAAAGACACGGTCTAGCAGACCGTCTATCGGCTAGGTCAAATAGGGTGCTTTGATATCAGCATGTCCAGCCTTAGAATTCAGTTCAGCGCGCTGGTCTGGGTCGAGATAAAATCACCAGTACCCAAGACCAGGCGGGCTCGCCGCGTTGGCTAATCCTGGTACATCTTGTAATCAATGTTCAGAAGAAAATCTGTGTTAGAGGGACGAGTCACCACGTACCAATAGCGACAACGATCGGTCGGACTATTCATCGTGGTGGTGACGCTCGGATTACGCGGGAAAGGTGCTTGTGTCCCGACAGGCTAGGATATAATGCTGAGGCGCTGCCCCAACCGTTCAGCGTGGGGTTTGCTACAACTTCCGAGTGCTACGTGTGCGAGACCATGTTATGTATGCACAAGGCCGACAATAGGACGTAGCCTTCGAGTTAGTACGTAGCGTGGTCGCACAAGCACAGTAGATCCTCCCCGCGCATCCTATTTATTAAGTTAATTCTATAGCAATACGATCACATGCGGATGGGCAGTGGCCGGTAGTCACACGCCTACCGCGGTGCTCAATGACCGGGACTAGAGAGGCGAAGATTATGGCGTGTGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGATTGCTTTCCCAATCTCCGAGCGATTTAGCGTGACAGCCCCAGGGAACCCACAAAATGCGATCGCAGTCCACCCGATCGTACACAGAAAGGAGGGTCCCCATACGCCGACGCACCTGTTCGCACGTCGTATGCATAAACGAGCCGCACGAACCAGAGAGCATAAAGAGGACCTCTAGCTCCTTTACAAAGTACAGGTTCGCCTGCCGCCGGGATGCCTTACCTAGACGCAATGACGGACGTATTCCTCTGGCCTCAACGGTTCCTGCTTCCGCTGGGATCCAAGATTGGCGGCCGAAGCCGCCTTTCCAAAGTGAGTCCTTCGTCTGTGACTAACTGTGCCAGATCGTCTTGCAAACTCCCGATCCAGTTTAACTCACCAAACTATAGCCGTACAGACCCAAATCTTAAGTCATATCACGCGACTAGCCTCTGCTCAATTTCTGTGCTCAAGGGTTTTGGTCCGCCCGAGCGGTGCAGCCGATTAGGACCATCTAATGCACTTGTTACAAGACTTCTTTTAAACACTTTCTTCCTGCCCAGTGGCGGATGATAATGGTTGTTGCCAGCCGGCGTGGAAGGTAACAGCACCGGTGCGAGCCTAATGTGCCGTCTCCACCAACACAGGGCTGTCCGGTCGTATAATAGGACTCCGCAATGGGGTTAGCAAGTGGCAGCCTAAACGATGTCGGGGACTCGCGATGTACATGCTCTGGTTCAATACATACGTGACCCGGCAGTTATCCTGCATCGGAACGTCAATCGTGCATCGGGCCAGCGTAATCGTGTCATCTGGGAGGCGGCCGTAGGATAAATAATTCAATAAAGATGTCGTTTTGCTAGTATACGCCTAGGCGTCACCCGCCATCTCTGTGCAGGTGGGCCGACGAGACACTGCCCCTGATTTCTCCGCTACTAATAGCACACACGGGGCAATACCAGCACAAGCCAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTGAGGCACGCCAATTCGCACTGTCGGGGTCGCTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTCCTTCCTTCCAGGGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCCCGTGCTTGCTCTCACGAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACCGGACAGCGCATCAACGCAACGGATCTACGTTACAGCGTGCATAGCGAGAACGGAGTTGCCGACGACGAAAGCGACACTGGGATCTGTCCGTCGTCATTCGCGGAAAGCATCCGCTCACGAGGCGGACACTGATTGACACGGTTTTGCAGAAGGTTAGGGGAATAGGTCAAATTGAGTGGCTTAAAAACGCTATGTCTGGGATTAAAGTGTAGTAAACTGCGGTCAACGGAGACGGTTTTAAGACAGGAGTTCGCAAAACCAGGCGGGGTCGCCACGACGGCTATTCCTGGTGGTTTAGGCGTACAATGTCCTGAAGAATATTTAAGAAAGAAGCACCCCTCGTCGCCTAGAATTACCTACCGCGGTCGACCATACCTTCGATTGTCGCGCCCACCCTCCCATTAGTCGGCAGAGGTGGTTGTGTTGCGATAGCCCAGCATGATATCCTAAGGCGTTACGCCGATGGATATCCCACGGAATTGCCATAGGCGCTGAACGCTACACGGACGATACGAACTTATGTATGGAGCGGGTCATCGAAAGGTCATACCCTTGTAGTTAACATGTAGCCCGGCCCTATTAGTACAGCAGTGCCTTGAGCGGCATTCTCATTATTAAGTTTTCTCTACAGCCAAACGACCAAGTGCACTTCCGCGGAGCGCGGTGGAGACTCGTCCACCCGGCAGCTCTGTAATAGGGACTAAAAGAGTGATGATAATCATGAGTGCCGCGTTATGGTGGTGTCGGAACAGAGCGGTCTTACGGCCAGTCGTATCCCTTCTCGAGTTCCGTCCGGTTAAGCGTGACACTCCCAGTGTACCCGCAAACCGTGATGGCTGTGCTTGGGGTCAATCGCATGTAGGATGGTCTCCAGACACCGGGGCACCAGTTTTCACGCCCAAAGCATAAACGACGAGCAGTCATGAGAGTCTTAGAACTGGACGTGCCGTTTCTCTGCGAACAACACCTCGAGCTGTACCGTTGTTGCGCTGCCTAGATGCAGTGCCGCTCCTATCACATTTGCCTCGACGACTGCCGCCTTCGCTGTTTCCCTAGACACTCAACAGTAAGCGCCTTTTGTAGGCAGGGGCACCCCCTGTCAGTGGCTGCGCCAAAACGTCTTCGGATCCCCTTGTCCAATCAAACTGACCGAATTCTTTCATTTAAGACCCTAATATGACATCATTAGTGACTAAATGCCACTCCCAAAATTCTGCCCAGAAGCGTTTAAGTTCGCCCCACTAAAGTTGTCTAAAACGA";
        let b = b"CTAAAGTGGCGAAATTTATGGTGTGTGACCCGTTATGCTCCATTTCGGTCAGTGGGTCATTGCTAGTAGTCGATTGCATTGTCATTCTCCGAGTGATTTAGCGTGACAGCCGCAGGGAACCCATAAAATGTAATCGTAGTCCATCTGATCGTACTTAGAAATGAAGGTCCCCTTTTACCCACGCACCTGTTTACTCGTCGTTTGCTTTTAAGAACCGCACGAACCACAGAGCATAAAGAGAACCTCTAGTTCCTTTACAAAGTACTGGTTCCCTTTTCAGCAAGATGCCTTATCTAAATGCAATGACAGACGTATTCCTCAGGCCACATCGCTTCCTACTTTCGCTGGGATCCATCATTGGCAGCTGAAACCGCCATTCCATAGTGAGTCCTTCGTCTGTGTCTTTCTGTGCCAAATCGTCTAGCAAATTGCTGATCCAGTTTATCTCACGAAATTATAGTCATACAGACCGAAATTTTAAATCAAATCACGCGACTAGGCTCAGCTTTATTTTAGTGGTCATGGGTTTTGGTCCGCCCGAGCGGTGCAACCGATTAGGACCATGTAAAACATTTGTTACAAGTCTTCTTTTAAATACAATCTTCCTGCTCAGTAGCGCATGATTATCGTTGTTGCTAGCCAGTGTGGTAAGTAACAGCACCACTGCGAGCCTAATGTGCCCTTTCCACGAACACAAGGCTATCCGATCCTATATTAGGATTCCGCAATGGGGTTAGCAAATCGCACCCTAAACGATATTGAAGACTTGCGATGTACATGCTTTGGTACAATACATACGTGTTCCAGTTGTTATCCTGTATCGGAACTTCAATTATGCATCGCACCAGCATATTCATGTCATCTAGGAAGAGCGCGTAGGATAAATAATTCAATTAAGATGTCGTTATGCTAGTATACGTCTACCCGTCACCGGCCATCTGTGTGCAGATGGGGCGACGAGTTATTGACCCTGATTTCTCCACTTCTAATACCACACACTGGGCAATACGAGCTCAAGCTAGTCTCGCAGTAACGCTCATCAGCTAACGAAAGAGTTAAAGGCTCGCTAATTCGCACTGTCAGGGTCTCTTGGGTGTTTTGCACTAGCGTCAGGTAGGCTAGTATGTGTTTTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAATGTGCAGCATACGTATTTGCTCGACGTGTTTAGTCTCTCATACTTCTCCTGGAGATCAAGGAAATGTTTCTTGTCCAAGTGGACAACGGTTCTACGGAATGGATCTACGTTACTGCCTGCATAAAGAAAACGGAGTTGCTAAGGACGAAAGCGACTTTAGGTTCTAACTGTTGACTTTGGCGGAAAAGTTTCATTCAGGAAGCAGACACTGATTGACACGGTTTAGCAGAACGTTTGAGGATTAGGTTAAATTGAGTGGTTTAATATTGGTATGTCTGGGATTAAAATATAGTATAGTGTGTTAATCGGAGACGAATTAAAGACACGAGTTCCCAAAATCAAGCGGGCTCATTACAACGGTTAATCCTGGTAGTTTACGTGAACAATGTTCTGAAGAAAATTTATGAAAAAAGGACCCGTCATCGCCTACAATTACCTACAACGGTCGACCATACCTTCGATTATCGTGGCCACTCTCGGATTACACGGCAGAGGTGGTTGTGTTCCGATAGGCCAGTATATTATTCTAAGGCGTTACCCTAATCATTTTTCATCGGATTTGCTATAGCCCTTGAACGCTACATGCACGAAACCAAATTATGTATACACTGGGTCATCAATAGGATATAGTCTTGTAGTTAACATGTAGCCCGGCCGTATTAGTACAGTAGAGCCTTCATTGACATTCTGTTTATTAAATTATTTCTACAGCAAAACGATCATATGCAAATCCACAGTGCGCGATAGAGATACATTCACTCGGCTGCTCTGTAATAGGGACTAAAAAAGTGATGATTATCATGAGTGCCCCGTTATGGTCGTGTTCGATCAGAGCGCTCTTACGAGCAGTCGTATACTTTCTCGAATTCCGTGCAGTTAAGCGTGACAGTCCCAGTGAACCCACAAAACGTGATGGCAGTCCATGCAATCATACGCAAGAAGGATGGTCTCCAGACACCGGCGCACCAGTTTTCACGCCGAAAGCATAAACGAGGAGCACAAATGAAAGTGTTTGAACTGGACCTGTAGTTTCTCTACGAAAAATACCTTGAGCTGTTGCGTTGTTGCGCTGCCTAGATGCAGTGTTGCACATATCACTTTTGCTTCAACGACTGCTGCTTTCGCTGTAACCCTAGACAGACAACAATAAGCGCTTTTTGTAGGCAAGAGCTCCGCCTATGACTAACTGCGCCAAAACATCTTCCAATCCCCTTATCCAATTTAATTCATCGAATTCTTACAATTTAGACCCTAATATCACATCATTAGACATTAATTGCCTCTGCCAAAATTCTGTCTACAAATGTTTTAGTTCGCTCCAGTAAAGTTGTTAATAACGACTACTAAATCCGCATGTTACGGGATTTCTTATTAATTCTTTTTTCGTAAGGAACAGCGGATCTTAATGGATGGCGCCAGGTGGTATGGAAGCTAATAGCGCGGGTGAGAGGGTAATTAGCCGTCTTCACCAACACAACGCTATCGGGTCATACTATAAGATTCCACAATGCGACTACTTATAAGATGTCTTAACGGTATCCGCAACTTGTGATGTGCCTACTATGCTTAAATGCATATCTCGCTCAGTAACTTTCCAATATGAGAGCATCAATTGTAGATCGGGCCGAGATAATCATGTCGTCACGGAACTTATTGTAAGAGTAATAATTTAAAAGAGATGTCAGTTTGCTGGTTCACGTAAAGGTTCCTCACACTACCTCTAAATAAGTGAGCGGTCGTGACATTATCCCTGATTTTCTCACTACTATTAGTACTCACGACACAATTCTACCACAGCCTTGTTTCGCCAGAATGCCAGTCAGCATAAAGAAGAGCTCAAGGCAGGTCAACTCGCATTGTGAGAGTTACATGAACGTTCGGCACTACCGACACGAACCTCAGTTAGCGTACATCCTACCAGAGGTCTGTGGCCCCGTGGTCAAAAGTGCGGATTTCGTATTTGCTGCTCGTCAGTACTTTCAGAATCATGACCTGCACGGTAAAAAGACGCTTATTATGGAGTTCGACATGGCAATAACGCGACGAATCTACGTCATGACGAGAATAGTATAAACAAAACTGCTGACGGCAGAAGCGTCAAAGAAGTCTGTAAATTGTTATTCGCGAAAAACATCCGTCTCCGTGGGGGATAATCACCGACGCCATTTTATAGAAGCCTAGGGGAACAGATTGGTTTAATTAGCTTAAGAAAGTAAATTCTGGGATTATACTGTAGTAATCACTAATTTACGGTGAGGGTTTTATGGCGGATTTTTACAAATTCAAACCAGGTGATTTCAACAAATTTTGTTGACGATTTAGGCGCACTATCCCCTAAACTACAAATTAAAAAATAGCGTTCCTTGACGGCTAGAATTACTTACCGGCCTTCACCATACCTTCGATATTCGCGCCCACTCTCCCATTAATCCGTACAAGTGGATGTAATGCGATTGTCCGCTAAGATATTCTAACGTGTAACGTAGATAAGTATTTTACAGAGTTGCCGTACGCGTTGAACACTTCACAGATGATAGGAATTTGCGTATAGAGCGTGTTATTGAGGAGTTATACACCCGTAGACTACAATGGGCCCAACTCAATCAGAACTCGAGTGCCTTGAATAACATACTCATCACTAAACATTCTCAACAATCAATCGAGCAAGTCCATTATCAACGAGTGTGTTGCAGTTTTATTCTCTTGCCAGCATTGTAATAGGCACTAAAAGAATGATGATAGTCATGAGTACTGAGCTAAGACGGCGTCGATGCATAGCGGACTTTCGGTCAATCACAATTCCTCACGAGACTCGTCCTGTTGAGCGTATCACTCTCAATGTACAAGCAACCCAAGAAGGCTGTGCCTGGACTCAACTGGATGCAGGATGAACTCCAGACACGGGGTCACTACTCTTCATACATAAAGCAAGAACGTCGAACAGTCATGAAAGTCTTAGTACCGCACGTACCATCTTACTGTGAATATTGCTTGAAGCTGTACCGTTATTGGGGGGCAAAGATGAAGTTCTCTTCTTTTCATAATTGTACTGACGACAGTCGTGTTCTCGGTTTCTTCAAAGGTTAAAGAATAAAGGCTTATTGTAGGCAGAGGAACGCCCTTTTAGTGGCTGGCGTTAAGTATCTTCGGACCCCCTTGTCTATCCAGATTAATCGAATTCTCTCATTTAGGACCTTAGTAAGTCATCATTGGTATTTGAATGCGACCTTGAAGAAACCGCTTAAAAATGTCAATGGTTGATCCACTAAACTTCATTTAATTAACTCCTAAATCAGCGCGATAGGCTATTAGAGGTTTAATTTTGTATAGCAAGGTACTTCCGATCTTAATGAATGGCCGGAAAAGGTACGGACGCGATATGCGAGGGTGAAAGGGCAAATAGACAGGTTCGTCTTTGTCACGCTAGGAGGCAATTCTATAAGAATGCATATTGCATCGATACATAAAATGTCTCGATCGCATGCGCAATTTGTGAAGTGTCTATTATCCCTAAGCCCATTTCCCGCATAATAACCCCTGATTGTATCCGCATTTGATGCTACCCAGGTTGAGTTAGCGTCGAGCTCGCGGAACTTATTGCATGAGTAGAGTTGAGTAAGAGCTGTTAGATGGCTCGCTGAACTAATAGTTGTCCACAGAACGTCAAGATTAGAAAACGGTTGTAGCATTATCGGAGGTTCTCTAACTACTATCAATACCCGTGTCTTGACTCTGCTGCGGCTACCTATCGCCTGAAAACCAGTTGGTGTTAAGGGATGCTCTGTCCAGGACGCCACATGTAGTGAAACTTACATGTTCGTTGGGTTCACCCGACT";
        // let a = b"ATCTAACTATTCCCTGTGCC";
        // let b = b"A";
        let score = |a: u8, b: u8| if a == b { 1 } else { -1 };
        let aligner = super::FastLSAAligner::with(a, b, score).with_block_size(150);
        let alignment = aligner.align();
        println!("FINAL ALIGNMENT:\n{alignment}");
    }
}
