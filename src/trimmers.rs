use std::collections::VecDeque;
use std::fmt;
use colored::Colorize;
use log::debug;
use crate::{color_qual_proportion, FastqRecord};
use crate::primers::PrimerSetMatch;

#[derive(PartialEq, Eq, PartialOrd, Debug, Clone, Ord)]
pub struct ReadSegment {
    pub start: usize,
    pub end: usize,
}

impl ReadSegment {
    pub fn join(segments: &Vec<ReadSegment>, intersect: &bool) -> Vec<ReadSegment> {
        if segments.len() == 0 {
            return Vec::new();
        }
        let mut segments = segments.clone();
        segments.sort();

        let mut accumulated_segments: Vec<ReadSegment> = vec![segments.get(0).unwrap().clone()];
        debug!("first accumulated_segments: {:?}", accumulated_segments);

        for segment in segments[1..segments.len()].iter() {
            debug!("accumulated_segments: {:?}, segment: {:?}", accumulated_segments, segment);
            let mut new_segments: Vec<ReadSegment> = Vec::new();
            for accumulated_segment in accumulated_segments.iter() {
                let set = match intersect {
                    &true => {ReadSegment::intersect_pair(accumulated_segment, segment)}
                    &false => { ReadSegment::join_pair(accumulated_segment, segment) }
                };
                match set {
                    Some(new_segment) => {
                        new_segments.push(new_segment);
                        debug!("new segments 1: {:?}", new_segments);
                    }
                    None => {
                        new_segments.push(accumulated_segment.clone());
                        new_segments.push(segment.clone());
                        debug!("new segments 2: {:?}", new_segments);
                    }
                }
            }
            new_segments.sort();
            accumulated_segments = new_segments;
        }
        accumulated_segments
    }


    pub fn join_pair(one: &ReadSegment, two: &ReadSegment) -> Option<ReadSegment> {
        let (leftmost, rightmost) = if std::cmp::min(one.start, two.start) == one.start { (one, two) } else { (two, one) };
        debug!("leftmost: {:?}, rightmost: {:?}", leftmost, rightmost);
        match (leftmost.end + 1 >= rightmost.start, leftmost.end >= rightmost.end) { // plus one for merging adjacent segments
            (true, true) => { Some(ReadSegment::new(leftmost.start, leftmost.end)) }
            (true, false) => { Some(ReadSegment::new(leftmost.start, rightmost.end)) }
            (false, true) => { panic!("impossible situation with two read segments") }
            (false, false) => { None }
        }
    }

    pub fn intersect_pair(one: &ReadSegment, two: &ReadSegment) -> Option<ReadSegment> {
        let (leftmost, rightmost) = if std::cmp::min(one.start, two.start) == one.start { (one, two) } else { (two, one) };
        debug!("leftmost: {:?}, rightmost: {:?}", leftmost, rightmost);
        match (leftmost.end + 1 >= rightmost.start, leftmost.end >= rightmost.end) { // plus one for merging adjacent segments
            (true, true) => { Some(ReadSegment::new(rightmost.start, rightmost.end)) }
            (true, false) => { Some(ReadSegment::new(rightmost.start, leftmost.end)) }
            (false, true) => { panic!("impossible situation with two read segments") }
            (false, false) => { None }
        }
    }

    pub fn contained(&self, index: &usize) -> bool {
        self.start <= *index && self.end > *index
    }
}

impl ReadSegment {
    pub fn new(start: usize, end: usize) -> ReadSegment {
        assert!(start <= end, "Read segment start must be less than or equal to end");
        ReadSegment { start, end }
    }
}

#[derive(Clone)]
pub struct TrimResult {
    keep_read: bool,
    remaining_read_segments: Vec<ReadSegment>,
}

impl std::fmt::Debug for TrimResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Keep: {}, segments {:?}", self.keep_read, self.remaining_read_segments)
    }
}

impl PartialEq for TrimResult {
    fn eq(&self, other: &Self) -> bool {
        self.keep_read == other.keep_read && self.remaining_read_segments.len() == other.remaining_read_segments.len() &&
            self.remaining_read_segments.iter().zip(other.remaining_read_segments.iter()).all(|(x, y)| x == y)
    }
}

impl TrimResult {
    pub fn keep(&self) -> bool {self.keep_read}

    pub fn from_read_segment(keep_read: bool, read_segments: ReadSegment) -> TrimResult {
        TrimResult { keep_read, remaining_read_segments: vec![read_segments] }
    }

    pub fn from_read(read: &FastqRecord) -> TrimResult {
        TrimResult::from_read_segment(true, ReadSegment::new(0, read.seq.len()))
    }

    pub fn join(trim_results: Vec<TrimResult>, interesct: &bool) -> TrimResult {
        let keep_read = trim_results.iter().map(|x| x.keep_read).all(|x| x);
        let read_segments = ReadSegment::join(&trim_results.iter().map(|x| x.remaining_read_segments.clone()).flatten().collect::<Vec<ReadSegment>>(), interesct);

        TrimResult { keep_read, remaining_read_segments: read_segments }
    }

    pub fn trim_results_to_reads(&self, read: &FastqRecord) -> Vec<FastqRecord> {
        let mut reads = Vec::new();
        if !self.keep_read {
            return reads;
        }

        for (index, read_segment) in self.remaining_read_segments.iter().enumerate() {
            let seq = read.seq[read_segment.start..read_segment.end].iter().map(|x| *x).collect();
            let qual = read.quals[read_segment.start..read_segment.end].iter().map(|x| *x).collect();

            match index {
                _x if _x > 0 => {
                    let new_name = format!("{}_segment{}", String::from_utf8(read.name.clone()).unwrap(), index);
                    reads.push(FastqRecord::new(new_name.into_bytes(), seq, qual));
                }
                _ => { reads.push(FastqRecord::new(read.name.clone(), seq.into(), qual.into())) }
            }
        }
        reads
    }


    pub fn print_format_read(&self, read: &FastqRecord) -> String {
        let mut seq_string = String::new();

        let mut current_slice_iter = self.remaining_read_segments.iter();
        let mut current_slice = current_slice_iter.next();

        for index in 0..read.seq.len() {
            let base = read.seq.get(index.clone()).unwrap();
            let quality = u32::from(read.quals[index]);
            assert!(quality >= 33, "Quality score must be at least 33. {} ", quality);

            let mut formated_base = format!("{}", char::from(*base));
            match *base {
                b'A' => { formated_base = formated_base.truecolor(0, color_qual_proportion(&quality), 0).to_string() }
                b'C' => { formated_base = formated_base.truecolor(color_qual_proportion(&quality), 0, 0).to_string() }
                b'G' => { formated_base = formated_base.truecolor(0, 0, color_qual_proportion(&quality)).to_string() }
                b'T' => { formated_base = formated_base.truecolor(color_qual_proportion(&quality), color_qual_proportion(&quality), 0).to_string() }
                _ => {
                    let qual_adjusted = color_qual_proportion(&quality);
                    formated_base = formated_base.truecolor(qual_adjusted.clone(), qual_adjusted.clone(), qual_adjusted.clone()).to_string()
                }
            }


            match current_slice {
                Some(slice) if slice.contained(&index) => {
                    formated_base = formated_base.on_truecolor(color_qual_proportion(&quality), color_qual_proportion(&quality), color_qual_proportion(&quality)).to_string();
                }
                _ => {
                    formated_base = formated_base.on_truecolor(255, 0, 0).to_string();
                }
            }

            if current_slice.is_some() && index > current_slice.unwrap().end {
                current_slice = current_slice_iter.next();
            }

            seq_string.push_str(&format!("{}", formated_base).to_string());
        }
        seq_string
    }
}

pub trait FastqTrimmer {
    fn trim(&self, read: &FastqRecord) -> TrimResult;
}


pub struct PrimerTrimmer {
    pub primer_set_matcher: PrimerSetMatch,
    // this describes where a primer is: if it's in one proportion worth at the front of the read,
    // then it's at the front, if it's in one proportion worth at the back of the read, then it's at the back,
    // otherwise it's in the middle and we currently kill the whole read
    pub end_proportion: f64,
    pub split_middle_adapters: bool,
}

impl PrimerTrimmer {
    pub fn new(primers: &Vec<Vec<u8>>, max_distance: &u8, end_proportion: &f64, split_middle_adapters: &bool) -> PrimerTrimmer {
        let mut all_orientations = primers.clone();
        all_orientations.extend(primers.iter().map(|x| PrimerTrimmer::rev_comp(x)));

        PrimerTrimmer {
            primer_set_matcher: PrimerSetMatch::new(&all_orientations, max_distance),
            end_proportion: end_proportion.clone(),
            split_middle_adapters: *split_middle_adapters,
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
        rev_comp.reverse();
        rev_comp
    }
}

impl FastqTrimmer for PrimerTrimmer {
    fn trim(&self, read: &FastqRecord) -> TrimResult {
        match self.primer_set_matcher.match_locations(&read.seq) {
            None => { TrimResult::from_read_segment(true, ReadSegment::new(0, read.seq.len())) }
            Some(x) => {
                let mut internal_clipped_primers = Vec::new();
                let front_start_pos = (self.end_proportion * read.seq.len() as f64).round() as usize;
                let end_start_pos = ((1.0 - self.end_proportion) * read.seq.len() as f64).round() as usize;
                let mut max_front = 0;
                let mut max_end = read.seq.len();

                for matching in x {
                    debug!("Matching {:?}", matching);
                    let match_start_end_pos = matching.1 + matching.0.len();

                    if match_start_end_pos < front_start_pos {
                        if match_start_end_pos > max_front {
                            debug!("max_front {:?}", match_start_end_pos);

                            max_front = match_start_end_pos;
                        }
                    } else if matching.1 >= end_start_pos {
                        if matching.1 < max_end {
                            debug!("max_end {:?}", matching.1);
                            max_end = matching.1;
                        }
                    } else {
                        debug!("MIDDLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        if self.split_middle_adapters {
                            internal_clipped_primers.push(ReadSegment::new(matching.1, matching.1 + matching.0.len()));
                        } else {
                            return TrimResult::from_read_segment(false, ReadSegment::new(0, read.seq.len()));
                        }
                    }
                }
                internal_clipped_primers.sort();

                if internal_clipped_primers.len() == 0 {
                    debug!("No clips {} {}",max_front, max_end);
                    return TrimResult::from_read_segment(true, ReadSegment::new(max_front, max_end));
                } else {
                    let mut resulting_read_segments = Vec::new();
                    internal_clipped_primers.iter().fold(max_front, |last_end, current| {
                        if current.start > last_end {
                            resulting_read_segments.push(ReadSegment::new(if last_end == 0 {0} else {last_end + 1}, current.start + 1));
                        }
                        current.end
                    });
                    resulting_read_segments.push(ReadSegment::new(internal_clipped_primers.get(internal_clipped_primers.len()-1).unwrap().end + 1, max_end));
                    debug!("Resulting read segments: {:?}", resulting_read_segments);
                    return TrimResult { keep_read: true, remaining_read_segments: resulting_read_segments };
                }
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
    fn trim(&self, read: &FastqRecord) -> TrimResult {
        TrimResult::from_read_segment(true, ReadSegment::new(0, BackTrimmer::trim_from_back(&read.quals, &self.window_min_qual_score, &self.window_size, &self.qual_score_base)))
    }
}

impl BackTrimmer {
    pub fn trim_from_back(quals: &Vec<u8>, window_min_qual_score: &u8, window_size: &u8, qual_score_base: &u8) -> usize {
        let range = (0..(quals.len() + 1 - *window_size as usize)).rev();
        let mut last_valid_trim = quals.len();
        for i in range {
            let average_qual = quals[i..(i + *window_size as usize)].iter().map(|&x| (x - *qual_score_base) as usize).sum::<usize>() as f64 / *window_size as f64;
            if average_qual >= f64::from(*window_min_qual_score) {
                return last_valid_trim;
            } else {
                last_valid_trim = i;
            }
        }
        last_valid_trim
    }
}

pub struct FrontBackTrimmer {
    pub window_size: u8,
    pub window_min_qual_score: u8,
    pub qual_score_base: u8,
}

impl FastqTrimmer for FrontBackTrimmer {
    fn trim(&self, read: &FastqRecord) -> TrimResult {
        TrimResult::from_read_segment(
            true,
            ReadSegment::new(read.quals.len() - BackTrimmer::trim_from_back(&read.quals.clone().into_iter().rev().collect(), &self.window_min_qual_score, &self.window_size, &self.qual_score_base),
                             BackTrimmer::trim_from_back(&read.quals, &self.window_min_qual_score, &self.window_size, &self.qual_score_base)))
    }
}

pub struct PolyXTrimmer {
    pub window_size: usize,
    pub minimum_g_proportion: f64,
    pub bases: Vec<u8>,
}

impl FastqTrimmer for PolyXTrimmer {
    fn trim(&self, read: &FastqRecord) -> TrimResult {
        let seq = &read.seq;
        let mut window: VecDeque<u8> = VecDeque::new();
        let mut found_g_window = false;
        debug!("looking for Gs");
        for index in (0..seq.len()).rev() {
            window.push_front(seq[index].clone());
            if window.len() > self.window_size {
                window.pop_back();
            }
            let g_prop = window.iter().map(|x| if self.bases.contains(x) { 1 } else { 0 }).sum::<usize>() as f64 / window.len() as f64;
            if g_prop < self.minimum_g_proportion && seq.len() - index > self.window_size {
                if found_g_window {
                    return TrimResult::from_read_segment(true, ReadSegment::new(0, index + 1));
                } else {
                    debug!("len()");
                    return TrimResult::from_read_segment(true, ReadSegment::new(0, seq.len()));
                }
            } else if g_prop >= self.minimum_g_proportion && window.len() == self.window_size {
                debug!("true");
                found_g_window = true;
            }
        };
        if found_g_window {
            TrimResult::from_read_segment(false, ReadSegment::new(0, 0))
        } else {
            TrimResult::from_read_segment(true, ReadSegment::new(0, seq.len()))
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::*;
    use crate::primers::PrimerSetMatch;
    use crate::trimmers::{PolyXTrimmer, ReadSegment, TrimResult};

    #[test]
    fn trim_two_color_off_the_end_test() {
        let g_trimmer = PolyXTrimmer {
            window_size: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b'G', b'g'],
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(false, ReadSegment::new(0, 0)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));


        let full_read_seq = "CGTACGCTAGACATTGTGCCGCATCATTAATGCAATGATTAGGAAATGACGGCTATACAAGGCAGATGAGGTTATTAGGGCATCCGCTTTAAGGCCGGTCCTACCAATGCGAACTCGATCATACCATGTGCACGCCTGTCTCTTAAACACATATGACGCTGCCGACGAGATCGTCCTCGTGTAGAACTCGGTGGACGCCGGATAATTAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".as_bytes().to_vec();
        let full_read_qual = (0..full_read_seq.len()).map(|_| b'C').into_iter().collect::<Vec<u8>>();
        let record = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: full_read_seq, quals: full_read_qual };
        assert_eq!(g_trimmer.trim(&record), TrimResult::from_read_segment(true, ReadSegment::new(0, 212)));

        let g_trimmer = PolyXTrimmer {
            window_size: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b'a', b'A'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 10)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGAAAAA".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));

        let g_trimmer = PolyXTrimmer {
            window_size: 5,
            minimum_g_proportion: 0.9,
            bases: vec![b't', b'T', b'a', b'A'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGG".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 10)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGATATA".as_bytes().to_vec(), quals: "CCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));

        let g_trimmer = PolyXTrimmer {
            window_size: 10,
            minimum_g_proportion: 0.9,
            bases: vec![b'g', b'G'],
        };
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GNCTTGGCACTCAGGCCAGGTTTGATCG".as_bytes().to_vec(), quals: "C#CCCCCCCCCCCCCCCCCCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 28)));
    }

    #[test]
    fn trim_the_back() {
        let g_trimmer = BackTrimmer {
            window_size: 5,
            window_min_qual_score: 5,
            qual_score_base: 32,
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "CCCCC!!!!!".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "!!!!!CCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 10)));


        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GNCTTGGCACTCAGGCCAGGTTTGATCG".as_bytes().to_vec(), quals: "C#CCCCCCCCCCCCCCCCCCCCCCCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 28)));
    }


    #[test]
    fn trim_adapters() {
        let g_trimmer = PrimerTrimmer {
            primer_set_matcher: PrimerSetMatch::new(&vec!["AAAAA".as_bytes().to_vec()], &1),
            end_proportion: 0.2,
            split_middle_adapters: false,
        };

        let fake_read = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAAGGGGGTTTTTCCCCCGGGGGTTTTTCCCCCGGGGGTTTTTCCCCC".as_bytes().to_vec(), quals: "AAAAAGGGGGTTTTTCCCCCGGGGGTTTTTCCCCCGGGGGTTTTTCCCCC".as_bytes().to_vec() };
        assert_eq!(g_trimmer.trim(&fake_read), TrimResult::from_read_segment(true, ReadSegment::new(5, 50)));

        let fake_read = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGTTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGTTTTTCCCCC".as_bytes().to_vec(), quals: "AAAAAGGGGGTTTTTCCCCCGGGGGTTTTTCCCCCGGGGGTTTTTCCCCC".as_bytes().to_vec() };
        assert_eq!(g_trimmer.trim(&fake_read), TrimResult::from_read_segment(false, ReadSegment::new(0, 50)));


        let g_trimmer = PrimerTrimmer {
            primer_set_matcher: PrimerSetMatch::new(&vec!["AAAAA".as_bytes().to_vec()], &1),
            end_proportion: 0.2,
            split_middle_adapters: true,
        };

        let fake_read = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGTTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGTTTTTCCCCC".as_bytes().to_vec(), quals: "AAAAAGGGGGTTTTTCCCCCGGGGGTTTTTCCCCCGGGGGTTTTTCCCCC".as_bytes().to_vec() };
        assert_eq!(g_trimmer.trim(&fake_read),
                   TrimResult{keep_read: true,
                       remaining_read_segments: vec![ReadSegment::new(0, 20), ReadSegment::new(25, 50)]});

        let fake_read = FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "AAAAACCCCCGGGGGCCCCCAAAAA".as_bytes().to_vec(), quals: "AAAAACCCCCGGGGGCCCCCAAAAA".as_bytes().to_vec() };
        let g_trimmer = PrimerTrimmer {
            primer_set_matcher: PrimerSetMatch::new(&vec!["GGGGG".as_bytes().to_vec()], &1),
            end_proportion: 0.1,
            split_middle_adapters: true,
        };
        assert_eq!(g_trimmer.trim(&fake_read),
                   TrimResult{keep_read: true,
                       remaining_read_segments: vec![ReadSegment::new(0, 10), ReadSegment::new(15, 25)]});
    }

    #[test]
    fn trim_front_and_back() {
        let g_trimmer = FrontBackTrimmer {
            window_size: 5,
            window_min_qual_score: 5,
            qual_score_base: 32,
        };

        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGG".as_bytes().to_vec(), quals: "CCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "CCCCC!!!!!".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(0, 5)));
        assert_eq!(g_trimmer.trim(&FastqRecord { name: "FAKEREAD".as_bytes().to_vec(), seq: "GGGGGGGGGG".as_bytes().to_vec(), quals: "!!!!!CCCCC".as_bytes().to_vec() }), TrimResult::from_read_segment(true, ReadSegment::new(5, 10)));
    }

    #[test]
    fn trim_result_join_check() {
        let trim_result1 = TrimResult::from_read_segment(true, ReadSegment::new(0, 5)); // keep the whole thing
        let trim_result2 = TrimResult::from_read_segment(true, ReadSegment::new(0, 5)); // keep the whole thing
        let trim_result3 = TrimResult::from_read_segment(true, ReadSegment::new(0, 20)); // keep the whole thing

        assert_eq!(TrimResult::join(vec![trim_result1.clone(), trim_result2.clone()], &false), trim_result1);
        assert_eq!(TrimResult::join(vec![trim_result1.clone(), trim_result2.clone(), trim_result3.clone()], &false), trim_result3);
    }

    #[test]
    fn join_pair() {
        let first = ReadSegment::new(0, 5);
        let second = ReadSegment::new(10, 50);
        let j = ReadSegment::join_pair(&first, &second);
        assert!(!j.is_some());

        let first = ReadSegment::new(0, 5);
        let second = ReadSegment::new(6, 50);
        let j = ReadSegment::join_pair(&first, &second);
        assert!(j.clone().is_some());
        assert_eq!(j.clone().unwrap().end, 50);
        assert_eq!(j.clone().unwrap().start, 0);

        let first = ReadSegment::new(6, 7);
        let second = ReadSegment::new(0, 50);
        let j = ReadSegment::join_pair(&first, &second);
        assert!(j.clone().is_some());
        assert_eq!(j.clone().unwrap().end, 50);
        assert_eq!(j.clone().unwrap().start, 0);

        let first = ReadSegment::new(6, 7);
        let second = ReadSegment::new(0, 6);
        let j = ReadSegment::join_pair(&first, &second);
        assert!(j.clone().is_some());
        assert_eq!(j.clone().unwrap().end, 7);
        assert_eq!(j.clone().unwrap().start, 0);

        let first = ReadSegment::new(31, 1495);
        let second = ReadSegment::new(0, 1495);
        let j = ReadSegment::intersect_pair(&first, &second);
        assert!(j.clone().is_some());
        assert_eq!(j.clone().unwrap().end, 1495);
        assert_eq!(j.clone().unwrap().start, 31);

    }

    #[test]
    fn join_vec() {
        let full_set = vec![ReadSegment::new(10, 50), ReadSegment::new(0, 5)];
        let j = ReadSegment::join(&full_set, &false);
        assert_eq!(j.len(), 2);
        assert_eq!(j[0].start, 0);
        assert_eq!(j[0].end, 5);
        assert_eq!(j[1].start, 10);
        assert_eq!(j[1].end, 50);

        let full_set = vec![ReadSegment::new(0, 5), ReadSegment::new(6, 50)];
        let j = ReadSegment::join(&full_set, &false);
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].start, 0);
        assert_eq!(j[0].end, 50);

        let full_set = vec![ReadSegment::new(0, 5), ReadSegment::new(4, 50)];
        let j = ReadSegment::join(&full_set, &false);
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].start, 0);
        assert_eq!(j[0].end, 50);

        let full_set = vec![ReadSegment::new(0, 50), ReadSegment::new(10, 20), ReadSegment::new(30, 40)];
        let j = ReadSegment::join(&full_set, &false);
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].start, 0);
        assert_eq!(j[0].end, 50);
    }

    #[test]
    fn read_segment_sort_check() {
        let first = ReadSegment::new(0, 5);
        let second = ReadSegment::new(10, 50);

        let mut vec = vec![second.clone(), first.clone()];
        vec.sort();
        assert_eq!(vec, vec![first, second]);

        let first = ReadSegment::new(6, 7);
        let second = ReadSegment::new(0, 50);

        let mut vec = vec![first.clone(), second.clone()];
        vec.sort();
        assert_eq!(vec, vec![second, first]);

        let first = ReadSegment::new(6, 7);
        let second = ReadSegment::new(0, 8);

        let mut vec = vec![first.clone(), second.clone()];
        vec.sort();
        assert_eq!(vec, vec![second, first]);
    }
}
