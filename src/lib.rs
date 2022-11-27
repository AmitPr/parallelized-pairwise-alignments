pub mod parallel;
pub mod sequential;
pub mod fast_lsa;

pub struct Alignment {
    pub score: i32,
    pub alignment: Vec<(char, char)>,
}
impl std::fmt::Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let s1 = self.alignment.iter().map(|(c, _)| c).collect::<String>();
        let s2 = self.alignment.iter().map(|(_, c)| c).collect::<String>();
        writeln!(f, "{s1}\n{s2}\nScore: {}", self.score)?;
        Ok(())
    }
}

pub trait Aligner<'a> {
    fn with(a: &'a [u8], b: &'a [u8], score: fn(u8, u8) -> i32) -> Self;
    fn align(&self) -> Alignment;
}

/// Return the last line of the Needleman-Wunsch matrix.
#[allow(clippy::needless_range_loop)]
pub fn nwscore(a: &[u8], b: &[u8], scoring: fn(u8, u8) -> i32) -> Vec<i32> {
    let mut score = (vec![0; b.len() + 1], vec![0; b.len() + 1]);
    for i in 1..=b.len() {
        score.0[i] = score.0[i - 1] + scoring(b'-', b[i - 1]);
    }
    for i in 0..a.len() {
        score.1[0] = score.0[0] + scoring(a[i], b'-');
        for j in 0..b.len() {
            let sub = score.0[j] + scoring(a[i], b[j]);
            let del = score.0[j + 1] + scoring(a[i], b'-');
            let ins = score.1[j] + scoring(b'-', b[j]);
            score.1[j + 1] = sub.max(del).max(ins);
        }
        score.0 = score.1.clone();
    }
    score.1
}
