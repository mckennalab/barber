use levenshtein_automata::{LevenshteinAutomatonBuilder, DFA};
use levenshtein_automata::SINK_STATE;

pub struct PrimerMatch {
    pub primer: Vec<u8>,
    pub lev_builder: LevenshteinAutomatonBuilder,
    pub dfa: DFA,
}

impl PrimerMatch {
    pub fn new(primer: &Vec<u8>, max_distance: u8) -> PrimerMatch {
        let lev_builder = LevenshteinAutomatonBuilder::new(max_distance, true);
        let dfa = lev_builder.build_dfa(&String::from_utf8(primer.clone()).unwrap().as_str());

        PrimerMatch {
            primer: primer.clone(),
            lev_builder,
            dfa,
        }
    }

    pub fn match_str(&self, seq: &Vec<u8>) -> bool {
        let mut state = self.dfa.initial_state();
        for &b in seq {
            state = self.dfa.transition(state, b);
            if state == SINK_STATE {
                return false;
            }
        }
        true
    }

    pub fn match_location(&self, seq: &Vec<u8>) -> Option<usize> {
        for i in 0..(seq.len() - self.primer.len()) {
            if self.match_str(&seq[i..i+self.primer.len()].to_vec()) {
                return Some(i.clone());
            }
        }
        None
    }
}

pub struct PrimerSetMatch {
    pub match_engines: Vec<PrimerMatch>,
}

impl PrimerSetMatch {
    pub fn new(primers: &Vec<Vec<u8>>, max_distance: &u8) -> PrimerSetMatch {
        let mut match_engines = Vec::new();
        for primer in primers {
            match_engines.push(PrimerMatch::new(primer, max_distance.clone()));
        }
        PrimerSetMatch {
            match_engines,
        }
    }

    pub fn match_locations(&self, seq: &Vec<u8>) -> Option<Vec<(&Vec<u8>,usize)>> {
        let mut matches = Vec::new();
        for match_engine in &self.match_engines{
            if let Some(location) = match_engine.match_location(seq) {
                matches.push((&match_engine.primer,location.clone()));
            }
        }
        if matches.len() > 0 {
            Some(matches)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::primers::{PrimerMatch, PrimerSetMatch};

    #[test]
    fn find_nanopore_primer() {
        let primer_seq = String::from("ACTTCGTTCAGTTACGTATTGCT");
        let partial_read = String::from("TGTCCTCTACTTCGTTCAGTTACGTATTGCTGTGAATTCGCCACCATGCTCGAGAGCTAGCGAATTCGAATTTGCCGGGTATGTGGCTGAAATCGAGTAACCGTTGCTAGGAGAGACCCTAGCAGAT");

        let matcher = PrimerMatch::new(&primer_seq.into_bytes(), 2);
        assert_eq!(matcher.match_location(&partial_read.into_bytes()), Some(6));

        let primer_seq = String::from("ACGTACGTACGTACGT");
        let partial_read = String::from("TGTCCTCTACTTCGTTCAGTTACGTATTGCTGTGAATTCGCCACCATGCTCGAGAGCTAGCGAATTCGAATTTGCCGGGTATGTGGCTGAAATCGAGTAACCGTTGCTAGGAGAGACCCTAGCAGAT");

        let matcher = PrimerMatch::new(&primer_seq.into_bytes(), 2);
        assert_eq!(matcher.match_location(&partial_read.into_bytes()), None);
    }

    #[test]
    fn any_primer_speed() {
        let primer_seq = String::from("ACTTCGTTCAGTTACGTATTGCT");
        let partial_read = String::from("TGTCCTCTACTTCGTTCAGTTACGTATTGCTGTGAATTCGCCACCATGCTCGAGAGCTAGCGAATTCGAATTTGCCGGGTATGTGGCTGAAATCGAGTAACCGTTGCTAGGAGAGACCCTAGCAGAT").into_bytes();

        let matcher = PrimerMatch::new(&primer_seq.into_bytes(), 2);
        for _j in 0..100000 {
            assert_eq!(matcher.match_location(&partial_read), Some(6));
        }
    }

    #[test]
    fn find_multiple_nanopore_primers() {
        let primer_seq = String::from("ACTTCGTTCAGTTACGTATTGCT");
        let primer_seq_nope = String::from("ACGTACGTACGTACGT");

        let partial_read = String::from("TGTCCTCTACTTCGTTCAGTTACGTATTGCTGTGAATTCGCCACCATGCTCGAGAGCTAGCGAATTCGAATTTGCCGGGTATGTGGCTGAAATCGAGTAACCGTTGCTAGGAGAGACCCTAGCAGAT");
        let not_a_good_partial_read = String::from("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        let matchers = PrimerSetMatch::new(&vec![primer_seq_nope.clone().into_bytes(), primer_seq.clone().into_bytes()], &2);
        assert_eq!(matchers.match_locations(&partial_read.into_bytes()), Some(vec![(&primer_seq.clone().into_bytes(),6)]));

        assert_eq!(matchers.match_locations(&not_a_good_partial_read.into_bytes()), None);
    }
}

