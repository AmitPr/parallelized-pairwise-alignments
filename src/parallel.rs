use std::sync::Arc;

use crate::{Aligner, Alignment};

// const PARALLEL_THRESHOLD: usize = 8;

pub struct ParallelHirschbergAligner<'a> {
    a: &'a [u8],
    b: &'a [u8],
    score: fn(u8, u8) -> i32,
}

impl<'a> Aligner<'a> for ParallelHirschbergAligner<'a> {
    fn with(a: &'a [u8], b: &'a [u8], score: fn(u8, u8) -> i32) -> Self {
        ParallelHirschbergAligner { a, b, score }
    }

    fn align(&self) -> Alignment {
        let score = self.score;
        if self.a.is_empty() {
            return Alignment {
                score: self.b.iter().map(|char| score(b'-', *char)).sum(),
                alignment: self.b.iter().map(|&b| ('-', b as char)).collect(),
            };
        }
        if self.b.is_empty() {
            return Alignment {
                score: self.a.iter().map(|char| score(*char, b'-')).sum(),
                alignment: self.a.iter().map(|&a| (a as char, '-')).collect(),
            };
        }
        if self.a.len() == 1 || self.b.len() == 1 {
            return self.align_nw();
        }
        let amid = self.a.len() / 2;
        let score_left = crate::nwscore(&self.a[..amid], self.b, score);
        let score_right = crate::nwscore(
            &self.a[amid..].iter().rev().cloned().collect::<Vec<_>>()[..],
            &self.b.iter().rev().cloned().collect::<Vec<_>>()[..],
            score,
        );

        let (ymid, _) = score_left
            .iter()
            .zip(score_right.iter().rev())
            .enumerate()
            .max_by_key(|(_, (l, r))| *l + *r)
            .unwrap();

        let a1 = Arc::new(&self.a[..amid]);
        let a2 = Arc::new(&self.a[amid..]);
        let b1 = Arc::new(&self.b[..ymid]);
        let b2 = Arc::new(&self.b[ymid..]);
        let score = Arc::new(score);

        // let mut aln1;
        // let mut aln2;
        // let do_aln1_parallel = a1.len() > PARALLEL_THRESHOLD && b1.len() > PARALLEL_THRESHOLD;
        // let do_aln2_parallel = a2.len() > PARALLEL_THRESHOLD && b2.len() > PARALLEL_THRESHOLD;
        // if  do_aln1_parallel && do_aln2_parallel {
        //     (aln1, aln2) = rayon::join(
        //         || ParallelHirschbergAligner::with(&a1, &b1, *score).align(),
        //         || ParallelHirschbergAligner::with(&a2, &b2, *score).align(),
        //     );
        // } else if do_aln1_parallel {
        //     (aln1, aln2) = rayon::join(
        //         || ParallelHirschbergAligner::with(&a1, &b1, *score).align(),
        //         || HirschbergAligner::with(&a2, &b2, *score).align(),
        //     );
        // }
        // else if do_aln2_parallel {
        //     (aln1, aln2) = rayon::join(
        //         || HirschbergAligner::with(&a1, &b1, *score).align(),
        //         || ParallelHirschbergAligner::with(&a2, &b2, *score).align(),
        //     );
        // }
        // else {
        //     (aln1, aln2) = rayon::join(
        //         || HirschbergAligner::with(&a1, &b1, *score).align(),
        //         || HirschbergAligner::with(&a2, &b2, *score).align(),
        //     );
        // }

        let (mut aln1, mut aln2) = rayon::join(
            || ParallelHirschbergAligner::with(&a1, &b1, *score).align(),
            || ParallelHirschbergAligner::with(&a2, &b2, *score).align(),
        );
        aln1.alignment.append(&mut aln2.alignment);
        aln1.score += aln2.score;
        aln1
    }
}

impl<'a> ParallelHirschbergAligner<'a> {
    /// Needleman-Wunsch alignment with one of the sequences being of length 1
    fn align_nw(&self) -> Alignment {
        let mut score = vec![vec![0; self.b.len() + 1]; self.a.len() + 1];
        let mut backtrack = vec![vec![0; self.b.len() + 1]; self.a.len() + 1];
        let scorer = self.score;
        for i in 1..=self.a.len() {
            score[i][0] = score[i - 1][0] + scorer(self.a[i - 1], b'-');
            backtrack[i][0] = 0;
        }
        for j in 1..=self.b.len() {
            score[0][j] = score[0][j - 1] + scorer(b'-', self.b[j - 1]);
            backtrack[0][j] = 1;
        }
        for i in 1..=self.a.len() {
            for j in 1..=self.b.len() {
                let scores = vec![
                    score[i - 1][j] + scorer(self.a[i - 1], b'-'), // deletion
                    score[i][j - 1] + scorer(b'-', self.b[j - 1]), // insertion
                    score[i - 1][j - 1] + scorer(self.a[i - 1], self.b[j - 1]), // match/mismatch
                ];
                let max_score = *scores.iter().max().unwrap();
                backtrack[i][j] = scores
                    .iter()
                    .enumerate()
                    .filter(|(_, &s)| s == max_score)
                    .map(|(i, _)| i)
                    .next()
                    .unwrap();
                score[i][j] = max_score;
            }
        }
        let mut i = self.a.len();
        let mut j = self.b.len();
        let mut alignment = Vec::new();
        while i > 0 || j > 0 {
            match backtrack[i][j] {
                0 => {
                    // alignment.insert(0, (self.a[i - 1] as char, '-'));
                    alignment.push((self.a[i - 1] as char, '-'));
                    i -= 1;
                }
                1 => {
                    alignment.push(('-', self.b[j - 1] as char));
                    // alignment.insert(0, ('-', self.b[j - 1] as char));
                    j -= 1;
                }
                2 => {
                    alignment.push((self.a[i - 1] as char, self.b[j - 1] as char));
                    // alignment.insert(0, (self.a[i - 1] as char, self.b[j - 1] as char));
                    i -= 1;
                    j -= 1;
                }
                _ => unreachable!(),
            }
        }
        alignment.reverse();
        Alignment {
            score: score[self.a.len()][self.b.len()],
            alignment,
        }
    }
}
