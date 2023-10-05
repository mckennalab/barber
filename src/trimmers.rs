use std::collections::VecDeque;
use crate::FastqRecord;

pub trait FastqTrimmer {
    fn trim(&self, read: &FastqRecord) -> (usize, usize);
}

pub struct BackTrimmer {
    pub window_size: u8,
    pub window_min_qual_score: u8,
    pub qual_score_base: u8,
}

impl FastqTrimmer for BackTrimmer {
    fn trim(&self, read: &FastqRecord) -> (usize, usize) {
        (0, BackTrimmer::trim_from_back(&read.quals, &self.window_min_qual_score, &self.window_size, &self.qual_score_base))
    }
}
impl BackTrimmer {
    pub fn trim_from_back(quals: &Vec<u8>, window_min_qual_score: &u8, window_size: &u8, qual_score_base: &u8) -> usize {
        let range = (0..(quals.len() + 1 - *window_size as usize)).rev();
        for i in range {
            let average_qual = quals[i..(i + *window_size as usize)].iter().map(|&x| (x - *qual_score_base) as usize).sum::<usize>() as f64 / *window_size as f64;
            if average_qual >= f64::from(*window_min_qual_score) {
                return i.clone() + (*window_size as usize)
            }
        }
        0
    }
}

pub struct FrontBackTrimmer {
    pub window_size: u8,
    pub window_min_qual_score: u8,
    pub qual_score_base: u8,
}

impl FastqTrimmer for FrontBackTrimmer {
    fn trim(&self, read: &FastqRecord) -> (usize, usize) {
        (read.quals.len() - BackTrimmer::trim_from_back(&read.quals.clone().into_iter().rev().collect(),&self.window_min_qual_score, &self.window_size, &self.qual_score_base),
         BackTrimmer::trim_from_back(&read.quals, &self.window_min_qual_score, &self.window_size, &self.qual_score_base))
    }
}

pub struct PolyXTrimmer {
    pub minimum_length: usize,
    pub minimum_g_proportion: f64,
    pub bases: Vec<u8>,
}

impl FastqTrimmer for PolyXTrimmer {
    fn trim(&self, read: &FastqRecord) -> (usize, usize) {
        let seq = &read.seq;
        let mut window: VecDeque<u8> = VecDeque::new();
        let mut found_g_window = false;
        //println!("looking for Gs");
        for index in (0..seq.len()).rev() {
            window.push_front(seq[index].clone());
            if window.len() > self.minimum_length {
                window.pop_back();
            }
            let g_prop = window.iter().map(|x| if self.bases.contains(x) { 1 } else { 0 }).sum::<usize>() as f64 / window.len() as f64;
            //println!("g_prop {} {} {}",g_prop,g_count as f64, (seq.len() - index));
            if g_prop < self.minimum_g_proportion && seq.len() - index > self.minimum_length {
                if found_g_window {
                    //println!("+1 {}",index + 1);
                    return (0, index + 1); // the last point that we were good
                } else {
                    //println!("len()");
                    return (0, seq.len()); // we never found a good window
                }
            } else if g_prop >= self.minimum_g_proportion {
                //println!("true");
                found_g_window = true;
            }
        };
        if found_g_window {
            (0, 0)
        } else {
            (0, seq.len())
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::*;
    use crate::trimmers::PolyXTrimmer;

    #[test]
    fn trim_two_color_off_the_end_test() {
        let g_trimmer = PolyXTrimmer {
            minimum_length: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b'G', b'g'],
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), (0,0));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (0,5));


        let full_read_seq = "CGTACGCTAGACATTGTGCCGCATCATTAATGCAATGATTAGGAAATGACGGCTATACAAGGCAGATGAGGTTATTAGGGCATCCGCTTTAAGGCCGGTCCTACCAATGCGAACTCGATCATACCATGTGCACGCCTGTCTCTTAAACACATATGACGCTGCCGACGAGATCGTCCTCGTGTAGAACTCGGTGGACGCCGGATAATTAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".as_bytes().to_vec();
        let full_read_qual = (0..full_read_seq.len()).map(|_| b'C').into_iter().collect::<Vec<u8>>();
        let record = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: full_read_seq, quals: full_read_qual };
        assert_eq!(g_trimmer.trim(&record), (0,212));

        let g_trimmer = PolyXTrimmer {
            minimum_length: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b'a', b'A'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (0,10));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGAAAAA".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (0,5));

        let g_trimmer = PolyXTrimmer {
            minimum_length: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b't', b'T', b'a', b'A'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (0,10));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGATATA".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (0,5));

    }
    #[test]
    fn trim_the_back() {
        let g_trimmer = BackTrimmer {
            window_size: 5,
            window_min_qual_score: 30,
            qual_score_base: 32,
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), (0,5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "CCCCC!!!!!".as_bytes().to_vec() }), (0,5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "!!!!!CCCCC".as_bytes().to_vec() }), (0,10));
        
    }


    #[test]
    fn trim_front_and_back() {
        let g_trimmer = FrontBackTrimmer {
            window_size: 5,
            window_min_qual_score: 30,
            qual_score_base: 32,
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), (0,5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "CCCCC!!!!!".as_bytes().to_vec() }), (0,5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "!!!!!CCCCC".as_bytes().to_vec() }), (5,10));

    }
    
}
