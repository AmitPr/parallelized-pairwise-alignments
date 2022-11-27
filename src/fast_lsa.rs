use std::sync::{Arc, Mutex};

use crate::{Aligner, Alignment};

const DEFAULT_BLOCKSIZE: usize = 10;

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

        let num_cached_rows = self.b.len() / self.block_size + 2;
        let num_cached_cols = self.a.len() / self.block_size + 2;

        let mut cached_rows = vec![vec![0; self.a.len() + 1]; num_cached_rows];
        let mut cached_cols = vec![vec![0; self.b.len() + 1]; num_cached_cols];

        let visited_blocks = Arc::new(Mutex::new(vec![
            vec![false; num_cached_cols - 1];
            num_cached_rows - 1
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

        todo!("implement fast lsa backtrace")
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
                    if row + 1 < num_cached_rows && !visited_blocks.lock().unwrap()[row + 1][col] {
                        let cached_rows = std::slice::from_raw_parts_mut(
                            cached_rows_ptr as *mut Vec<i32>,
                            num_cached_rows - row - 1,
                        );
                        let cached_cols = std::slice::from_raw_parts_mut(
                            cached_cols_ptr as *mut Vec<i32>,
                            num_cached_cols - row,
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
                    if col + 1 < num_cached_cols && !visited_blocks_2.lock().unwrap()[row][col + 1]
                    {
                        let cached_rows = std::slice::from_raw_parts_mut(
                            cached_rows_ptr as *mut Vec<i32>,
                            num_cached_rows - row,
                        );
                        let cached_cols = std::slice::from_raw_parts_mut(
                            cached_cols_ptr as *mut Vec<i32>,
                            num_cached_cols - col - 1,
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
}

#[cfg(test)]
mod tests {
    use crate::Aligner;

    #[test]
    fn doit() {
        let a = b"CCGTTATGCTCCATTTCGG";
        let b = b"CCTGCGCTATCTGCCTGTC";
        let score = |a: u8, b: u8| if a == b { 1 } else { -1 };
        let aligner = super::FastLSAAligner::with(a, b, score);
        let alignment = aligner.align();
    }
}
