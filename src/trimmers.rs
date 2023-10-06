use std::collections::VecDeque;
use crate::FastqRecord;
use crate::primers::PrimerSetMatch;

pub trait FastqTrimmer {
    fn trim(&self, read: &FastqRecord) -> (bool, usize, usize);
}


pub struct PrimerTrimmer {
    pub primer_set_matcher: PrimerSetMatch,
    pub max_distance: u8,
    // this describes where a primer is: if it's in one proportion worth at the front of the read,
    // then it's at the front, if it's in one proportion worth at the back of the read, then it's at the back,
    // otherwise it's in the middle and we currently kill the whole read
    pub end_proportion: f64,
}

impl PrimerTrimmer {
    pub fn new(primers: &Vec<Vec<u8>>, max_distance: &u8, end_proportion: &f64) -> PrimerTrimmer {
        let mut all_orientations = primers.clone();
        all_orientations.extend(primers.iter().map(|x| PrimerTrimmer::rev_comp(x)));

        PrimerTrimmer {
            primer_set_matcher: PrimerSetMatch::new(&all_orientations, max_distance),
            max_distance: max_distance.clone(),
            end_proportion: end_proportion.clone(),
        }
    }

    pub fn rev_comp(seq: &Vec<u8>) -> Vec<u8> {
        let mut rev_comp = Vec::new();
        for base in seq.iter().rev() {
            rev_comp.push(match base {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                _ => panic!("Invalid base in primer"),
            });
        }
        rev_comp
    }
}

impl FastqTrimmer for PrimerTrimmer {
    fn trim(&self, read: &FastqRecord) -> (bool, usize, usize) {
        match self.primer_set_matcher.match_locations(&read.seq) {
            None => { (true, 0, read.seq.len()) }
            Some(x) => {
                let front_start_pos = (self.end_proportion * read.seq.len() as f64).round() as usize;
                let end_start_pos = ((1.0 - self.end_proportion) * read.seq.len() as f64).round() as usize;
                let mut max_front = 0;
                let mut max_end = read.seq.len();

                for matching in x {
                    if matching.1 < front_start_pos as usize {
                        if matching.1 > max_front {
                            max_front = matching.1;
                        }
                    } else if matching.1 >= end_start_pos as usize {
                        if matching.1 < max_end {
                            max_end = matching.1;
                        }
                    } else {
                        return (false, 0, read.seq.len());
                    }
                }
                (true, max_front, max_end)
            }
        }
    }
}

pub struct BackTrimmer {
    pub window_size: u8,
    pub window_min_qual_score: u8,
    pub qual_score_base: u8,
}

impl FastqTrimmer for BackTrimmer {
    fn trim(&self, read: &FastqRecord) -> (bool, usize, usize) {
        (true, 0, BackTrimmer::trim_from_back(&read.quals, &self.window_min_qual_score, &self.window_size, &self.qual_score_base))
    }
}

impl BackTrimmer {
    pub fn trim_from_back(quals: &Vec<u8>, window_min_qual_score: &u8, window_size: &u8, qual_score_base: &u8) -> usize {
        let range = (0..(quals.len() + 1 - *window_size as usize)).rev();
        for i in range {
            let average_qual = quals[i..(i + *window_size as usize)].iter().map(|&x| (x - *qual_score_base) as usize).sum::<usize>() as f64 / *window_size as f64;
            if average_qual >= f64::from(*window_min_qual_score) {
                return i.clone() + (*window_size as usize);
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
    fn trim(&self, read: &FastqRecord) -> (bool, usize, usize) {
        (true,
         read.quals.len() - BackTrimmer::trim_from_back(&read.quals.clone().into_iter().rev().collect(), &self.window_min_qual_score, &self.window_size, &self.qual_score_base),
         BackTrimmer::trim_from_back(&read.quals, &self.window_min_qual_score, &self.window_size, &self.qual_score_base))
    }
}

pub struct PolyXTrimmer {
    pub window_size: usize,
    pub minimum_g_proportion: f64,
    pub bases: Vec<u8>,
}

impl FastqTrimmer for PolyXTrimmer {
    fn trim(&self, read: &FastqRecord) -> (bool, usize, usize) {
        let seq = &read.seq;
        let mut window: VecDeque<u8> = VecDeque::new();
        let mut found_g_window = false;
        //println!("looking for Gs");
        for index in (0..seq.len()).rev() {
            window.push_front(seq[index].clone());
            if window.len() > self.window_size {
                window.pop_back();
            }
            let g_prop = window.iter().map(|x| if self.bases.contains(x) { 1 } else { 0 }).sum::<usize>() as f64 / window.len() as f64;
            if g_prop < self.minimum_g_proportion && seq.len() - index > self.window_size {
                if found_g_window {
                    return (true, 0, index + 1); // the last point that we were good
                } else {
                    //println!("len()");
                    return (true, 0, seq.len()); // we never found a good window
                }
            } else if g_prop >= self.minimum_g_proportion && window.len() == self.window_size {
                //println!("true");
                found_g_window = true;
            }
        };
        if found_g_window {
            (false, 0, 0)
        } else {
            (true, 0, seq.len())
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
            window_size: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b'G', b'g'],
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), (false, 0, 0));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 5));


        let full_read_seq = "CGTACGCTAGACATTGTGCCGCATCATTAATGCAATGATTAGGAAATGACGGCTATACAAGGCAGATGAGGTTATTAGGGCATCCGCTTTAAGGCCGGTCCTACCAATGCGAACTCGATCATACCATGTGCACGCCTGTCTCTTAAACACATATGACGCTGCCGACGAGATCGTCCTCGTGTAGAACTCGGTGGACGCCGGATAATTAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".as_bytes().to_vec();
        let full_read_qual = (0..full_read_seq.len()).map(|_| b'C').into_iter().collect::<Vec<u8>>();
        let record = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: full_read_seq, quals: full_read_qual };
        assert_eq!(g_trimmer.trim(&record), (true, 0, 212));

        let g_trimmer = PolyXTrimmer {
            window_size: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b'a', b'A'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 10));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGAAAAA".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 5));

        let g_trimmer = PolyXTrimmer {
            window_size: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b't', b'T', b'a', b'A'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 10));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGATATA".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 5));

        let g_trimmer = PolyXTrimmer {
            window_size: 10,
            minimum_g_proportion: 0.9,
            bases: vec![b'g', b'G'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GNCTTGGCACTCAGGCCAGGTTTGATCG".as_bytes().to_vec(), quals: "C#CCCCCCCCCCCCCCCCCCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 28));
    }

    #[test]
    fn trim_the_back() {
        let g_trimmer = BackTrimmer {
            window_size: 5,
            window_min_qual_score: 30,
            qual_score_base: 32,
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), (true, 0, 5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "CCCCC!!!!!".as_bytes().to_vec() }), (true, 0, 5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "!!!!!CCCCC".as_bytes().to_vec() }), (true, 0, 10));


        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GNCTTGGCACTCAGGCCAGGTTTGATCG".as_bytes().to_vec(), quals: "C#CCCCCCCCCCCCCCCCCCCCCCCCCC".as_bytes().to_vec() }), (true, 0, 28));
    }


    #[test]
    fn trim_front_and_back() {
        let g_trimmer = FrontBackTrimmer {
            window_size: 5,
            window_min_qual_score: 30,
            qual_score_base: 32,
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), (true, 0, 5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "CCCCC!!!!!".as_bytes().to_vec() }), (true, 0, 5));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "!!!!!CCCCC".as_bytes().to_vec() }), (true, 5, 10));
    }
}
