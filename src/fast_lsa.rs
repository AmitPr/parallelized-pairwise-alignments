use std::sync::{Arc, Mutex};

use crate::{Aligner, Alignment};

const DEFAULT_BLOCKSIZE: usize = 150;

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

        let mut cached_rows = vec![vec![0; self.a.len() + 1]; num_cached_rows];
        let mut cached_cols = vec![vec![0; self.b.len() + 1]; num_cached_cols];

        let visited_blocks = Arc::new(Mutex::new(vec![
            vec![false; num_cached_cols];
            num_cached_rows
        ]));

        /*
         * 1. Fill in the cache grid
         */

        // First row and first column are always just insertion and deletion scores
        self.a.iter().enumerate().for_each(|(i, char)| {
            cached_rows[0][i + 1] = cached_rows[0][i] + score(*char, b'-');
        });
        self.b.iter().enumerate().for_each(|(i, char)| {
            cached_cols[0][i + 1] = cached_cols[0][i] + score(b'-', *char);
        });

        self.forward(
            0,
            0,
            cached_rows.as_mut_slice(),
            cached_cols.as_mut_slice(),
            score,
            visited_blocks,
        );

        /*
         * 2. Backtrace
         */
        let mut last_a = &self.a[(num_grid_cols - 1) * self.block_size..];
        let mut last_b = &self.b[(num_grid_rows - 1) * self.block_size..];
        let mut start_row = &cached_rows[num_grid_rows][self.a.len() - last_a.len()..];
        let mut start_col = &cached_cols[num_grid_cols][self.b.len() - last_b.len()..];
        let mut start_idx = (last_a.len(), last_b.len());
        let mut next_row = num_grid_rows - 1;
        let mut next_col = num_grid_cols - 1;

        let mut alignment = vec![];
        let mut alignment_score = cached_rows.last().unwrap().last().unwrap();

        loop {
            let result = Self::backtrace_grid_cell(
                last_a, last_b, start_row, start_col, start_idx, next_row, next_col, score,
            );
            alignment.splice(0..0, result.0.alignment);

            if next_row == 0 && next_col == 0 {
                break;
            }
            // (sub_alignment, start_idx, next_row, next_col)
            if next_row > result.2 {
                start_idx = (self.block_size, result.1);
            } else {
                start_idx = (result.1, self.block_size);
            }
            next_row = result.2;
            next_col = result.3;

            last_a = &self.a
                [next_col * self.block_size..((next_col + 1) * self.block_size).min(self.a.len())];
            last_b = &self.b
                [next_row * self.block_size..((next_row + 1) * self.block_size).min(self.b.len())];

            start_row = &cached_rows[next_row][next_col * self.block_size
                ..((next_col + 1) * self.block_size + 1).min(self.a.len() + 1)];
            start_col = &cached_cols[next_col][next_row * self.block_size
                ..((next_row + 1) * self.block_size + 1).min(self.b.len() + 1)];
        }

        Alignment {
            alignment,
            score: *alignment_score,
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
        row: usize,
        col: usize,
        cached_rows: &mut [Vec<i32>],
        cached_cols: &mut [Vec<i32>],
        score: fn(u8, u8) -> i32,
        visited_blocks: Arc<Mutex<Vec<Vec<bool>>>>,
    ) {
        let num_grid_rows = (self.b.len() + self.block_size - 1) / self.block_size;
        let num_grid_cols = (self.a.len() + self.block_size - 1) / self.block_size;
        assert!(row < num_grid_rows, "row out of bounds");
        assert!(col < num_grid_cols, "col out of bounds");
        visited_blocks.lock().unwrap()[row][col] = true;

        let row_start = col * self.block_size;
        let row_end = ((col + 1) * self.block_size).min(self.b.len());
        let col_start = row * self.block_size;
        let col_end = ((row + 1) * self.block_size).min(self.a.len());

        // Correctness: We never mutably borrow the same parts of a row or column twice
        unsafe {
            let start_row = cached_rows.get_unchecked(row)[row_start..=row_end].as_ptr();
            let start_col = cached_cols.get_unchecked(col)[col_start..=col_end].as_ptr();
            let output_row =
                cached_rows.get_unchecked_mut(row + 1)[row_start..=row_end].as_mut_ptr();
            let output_col =
                cached_cols.get_unchecked_mut(col + 1)[col_start..=col_end].as_mut_ptr();
            // convert pointers to slices
            let start_row = std::slice::from_raw_parts(start_row, row_end - row_start + 1);
            let start_col = std::slice::from_raw_parts(start_col, col_end - col_start + 1);
            let output_row = std::slice::from_raw_parts_mut(output_row, row_end - row_start + 1);
            let output_col = std::slice::from_raw_parts_mut(output_col, col_end - col_start + 1);

            Self::fill_grid_block(
                &self.a[row_start..row_end],
                &self.b[col_start..col_end],
                start_row,
                start_col,
                output_row,
                output_col,
                score,
            );
            println!("{}x{}: Input row: {:?}, Output row: {:?}", row, col, start_row, output_row);
        }
        // Correctness: We can assure ourselves that cached_rows and cached_cols are not
        // dropped before the forward call is finished. We can also assure ourselves that
        // the mutable pointers that we pass are not used to mutate the same parts of the
        // cache vectors twice.
        unsafe {
            let num_cached_rows = cached_rows.len();
            let num_cached_cols = cached_cols.len();
            let cached_rows_ptr = cached_rows.as_mut_ptr() as usize;
            let cached_cols_ptr = cached_cols.as_mut_ptr() as usize;

            let visited_blocks_2 = Arc::clone(&visited_blocks);
            rayon::join(
                || {
                    if row + 1 < num_grid_rows
                        && !visited_blocks
                            .lock()
                            .unwrap_or_else(|_| panic!("Could not lock mutex"))[row + 1][col]
                    {
                        let cached_rows = std::slice::from_raw_parts_mut(
                            cached_rows_ptr as *mut Vec<i32>,
                            num_cached_rows,
                        );
                        let cached_cols = std::slice::from_raw_parts_mut(
                            cached_cols_ptr as *mut Vec<i32>,
                            num_cached_cols,
                        );
                        self.forward(
                            row + 1,
                            col,
                            cached_rows,
                            cached_cols,
                            score,
                            visited_blocks,
                        );
                    }
                },
                || {
                    if col + 1 < num_grid_cols
                        && !visited_blocks_2
                            .lock()
                            .unwrap_or_else(|_| panic!("Could not lock mutex"))[row][col + 1]
                    {
                        {
                            let cached_rows = std::slice::from_raw_parts_mut(
                                cached_rows_ptr as *mut Vec<i32>,
                                num_cached_rows,
                            );
                            let cached_cols = std::slice::from_raw_parts_mut(
                                cached_cols_ptr as *mut Vec<i32>,
                                num_cached_cols,
                            );
                            self.forward(
                                row,
                                col + 1,
                                cached_rows,
                                cached_cols,
                                score,
                                visited_blocks_2,
                            );
                        }
                    }
                },
            );
        }
    }

    pub fn fill_grid_block(
        a: &[u8],
        b: &[u8],
        start_row: &[i32],
        start_col: &[i32],
        output_row: &mut [i32],
        output_col: &mut [i32],
        score: fn(u8, u8) -> i32,
    ) {
        let n = a.len();
        let m = b.len();
        // first index already computed in start_{row,col}
        let mut scores = (vec![0; n + 1], vec![0; n + 1]);
        scores.0 = start_row.to_vec();
        output_col[0] = start_row[n];
        for i in 0..m {
            scores.1[0] = start_col[i + 1];

            for j in 0..n {
                scores.1[j + 1] = {
                    let subst = scores.0[j] + score(a[j], b[i]);
                    let del = scores.1[j] + score(a[j], b'-');
                    let ins = scores.0[j + 1] + score(b'-', b[i]);
                    subst.max(del).max(ins)
                };
                // Fill output column
                if j == n - 1 {
                    output_col[i + 1] = scores.1[j + 1];
                }
            }
            scores.0 = scores.1.clone();
        }
        // Fill output row
        output_row.copy_from_slice(&scores.1);
    }

    /// Returns: (alignment, next_start_idx, next_row, next_col)
    pub fn backtrace_grid_cell(
        a: &[u8],
        b: &[u8],
        start_row: &[i32],
        start_col: &[i32],
        start_idx: (usize, usize),
        row: usize,
        col: usize,
        score: fn(u8, u8) -> i32,
    ) -> (Alignment, usize, usize, usize) {
        // println!("Backward: row: {}, col: {}", row, col);
        let n = a.len();
        let m = b.len();
        let mut scores = vec![vec![0; n + 1]; m + 1];
        let mut backtrack = vec![vec![0; n + 1]; m + 1];
        if row == 0 {
            for j in 0..n {
                backtrack[0][j + 1] = 2;
            }
        }
        if col == 0 {
            for i in 0..m {
                backtrack[i + 1][0] = 1;
            }
        }
        scores[0] = start_row.to_vec();
        for i in 0..m {
            scores[i + 1][0] = start_col[i + 1];
            for j in 0..n {
                let nexts = vec![
                    scores[i][j] + score(a[j], b[i]),     // substitution
                    scores[i][j + 1] + score(b'-', b[i]), // deletion (go down)
                    scores[i + 1][j] + score(a[j], b'-'), // insertion (go right)
                ];
                let max_score = *nexts.iter().max().unwrap();
                backtrack[i + 1][j + 1] = nexts
                    .iter()
                    .enumerate()
                    .filter(|(_, &s)| s == max_score)
                    .map(|(i, _)| i)
                    .next()
                    .unwrap();
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
        while condition(i, j) {
            match backtrack[i][j] {
                0 => {
                    alignment.push((a[j - 1] as char, b[i - 1] as char));
                    i -= 1;
                    j -= 1;
                }
                1 => {
                    alignment.push(('-', b[i - 1] as char));
                    i -= 1;
                }
                2 => {
                    alignment.push((a[j - 1] as char, '-'));
                    j -= 1;
                }
                _ => unreachable!(),
            }
        }
        alignment.reverse();
        let alignment = Alignment {
            score: scores[start_idx.0][start_idx.1],
            alignment,
        };
        if row == 0 && col == 0 {
            (alignment, 0, 0, 0)
        } else {
            let next_start_idx = if i > 0 { i } else { j };
            let next_row = if i > 0 { row } else { row - 1 };
            let next_col = if i > 0 { col - 1 } else { col };
            (alignment, next_start_idx, next_row, next_col)
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
        let score = |a: u8, b: u8| if a == b { 1 } else { -1 };
        let aligner = super::FastLSAAligner::with(a, b, score);
        let alignment = aligner.align();
        println!("FINAL ALIGNMENT:\n{}", alignment);
    }
}
