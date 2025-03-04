extern crate ndarray;
//extern crate rustworkx;
extern crate bio;
extern crate numpy;
extern crate petgraph;
extern crate rustworkx_core;

use bio::alignment::pairwise;
use bio::alignment::AlignmentOperation as Op;
use bio::alphabets;
use core::cell::Cell;
use ndarray::prelude::*;
use numpy::PyArray2;
use petgraph::graph::Graph;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::pyfunction;
use pyo3::PyErr;
use pyo3::PyResult;
use pyo3::Python;
use rustworkx_core::connectivity::connected_components;
use std::cmp;
use std::collections::HashMap;
use std::collections::HashSet;

type Matrix = Array2<u16>;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

fn to_graph<'a, N, E>(nodes: &Vec<Vec<usize>>) -> Graph<N, E>
where
    N: Default,
    E: Default,
{
    let mut edges: Vec<(u32, u32)> = vec![];
    for l in nodes {
        for p in l.windows(2) {
            edges.push((p[0] as u32, p[1] as u32));
        }
    }
    Graph::from_edges(&edges)
}

#[pyfunction]
fn find_repeats(sseq: String) -> PyResult<HashMap<usize, Vec<Vec<usize>>>> {
    match rust_find_repeats(sseq) {
        None => Err(PyErr::new::<PyTypeError, _>("foo")),
        Some(m) => Ok(m),
    }
}

fn rust_find_repeats(sseq: String) -> Option<HashMap<usize, Vec<Vec<usize>>>> {
    let seq: Vec<u8> = sseq.bytes().collect();
    let repeat_lower_limit: u16 = 5;
    let threshold_str2bp: usize = 10;
    let threshold_str3bp: usize = 6;

    if seq.len() < threshold_str2bp * 2 + 2 {
        return None;
    }

    let mut m = Matrix::from_elem((seq.len(), seq.len()), 0);
    // fill first row of matrix
    let first = seq.iter().next().unwrap();
    for (i, b) in seq.iter().enumerate().skip(1) {
        if b == first {
            m[[0, i]] = 1;
        }
    }
    // fill the rest of the diagonal
    for i in 1..seq.len() {
        for j in i + 1..seq.len() {
            if seq[i] == seq[j] {
                m[[i, j]] = m[[i - 1, j - 1]] + 1;
            }
        }
    }

    let mut d: HashMap<usize, Vec<Vec<usize>>> = HashMap::new();
    d.insert(2, vec![]);
    // find 2bp STR, (ab)n pattern
    for i in 1..seq.len() - (2 * threshold_str2bp) + 2 {
        // this is a self repeat, so second instance will always begin 2 units from diagonal
        if m[[i, i + 2]] == 2 {
            let j = i + 2;
            let mut a = 2;
            let mut b = 0;
            let mut n = 1;
            // require a=2 to indicate local repeat 2bp long; require b=0 to prevent counting homopolymers
            while b == 0 && a == 2 && (j + n * 2) < seq.len() {
                a = m[[i, j + n * 2]];
                b = m[[i, j + n * 2 - 1]];
                if a == 2 && b == 0 {
                    n += 1;
                }
            }
            // if above limit, record in array as [pos (0-based), number of repeats]
            if n >= (threshold_str2bp - 1) {
                let mut str_tuple = vec![(i - 1) as usize];
                for x in 1..n + 1 {
                    str_tuple.push((i - 1 + x * 2) as usize);
                }
                d.entry(2).or_insert_with(|| vec![]).push(str_tuple);
            }
        }
    }

    // find 3bp STR, (aab)n, (abb)n, or (abc)n pattern
    d.insert(3, vec![]);
    for i in 1..seq.len() - (3 * threshold_str3bp) {
        // this is a self repeat, so second instance will always begin 3 units from diagonal
        if m[[i, i + 3]] == 3 {
            let j = i + 3;
            let mut a = 3;
            let mut b = 0;
            let mut n = 1;
            // require a=3 to indicate local repeat 3bp long; require b=0 at position 2,3 of repeat to prevent counting homopolymers
            while b == 0 && a == 3 && (j + n * 3) < seq.len() {
                a = m[[i, j + n * 3]];
                b = cmp::min(m[[i, j + n * 3 - 1]], m[[i, j + n * 3 - 2]]);
                if a == 3 && b == 0 {
                    n += 1;
                }
            }
            // if above limit, record in array as [pos (0-based), number of repeats]
            if n >= (threshold_str3bp - 1) {
                let mut str_tuple = vec![(i - 2) as usize];
                for x in 1..n + 1 {
                    str_tuple.push((i - 2 + x * 3) as usize);
                }
                d.entry(3).or_insert_with(|| vec![]).push(str_tuple);
            }
        }
    }

    for i in 1..seq.len() {
        // don't consider lower half of array, diagonal (self repeat), or 1 base offset (homopolymer repeat)
        for j in i + 2..seq.len() {
            // don't count subunits of a repeat that is already counted
            if m[[i, j]] >= repeat_lower_limit
                && m[[i, j]] > *m.slice(s![i, ..j]).iter().max().unwrap()
            {
                // not an edge, a default repeat end
                if i < seq.len() - 1 && j < seq.len() - 1 {
                    if m[[i, j]] > m[[i + 1, j + 1]]
                        && m[[i, j]] != m[[i - 2, j]] + 2
                        && m[[i, j]] != m[[i - 3, j]] + 3
                    {
                        let rep_len = m[[i, j]] as usize;
                        let mut rep_pos =
                            vec![(i - rep_len + 1) as usize, (j - rep_len + 1) as usize];
                        // look for other instances of the same repeat, not str2 or str3
                        for k in i + 4..seq.len() {
                            if m[[i, k]] == m[[i, j]] {
                                rep_pos.push((k - rep_len + 1) as usize)
                            }
                        }
                        d.entry(rep_len).or_insert_with(|| vec![]).push(rep_pos);
                    }
                // add repeats at sequence edges
                } else {
                    if m[[i, j]] != m[[i - 2, j]] + 2 && m[[i, j]] != m[[i - 3, j]] + 3 {
                        let rep_len = m[[i, j]] as usize;
                        let mut rep_pos =
                            vec![(i - rep_len + 1) as usize, (j - rep_len + 1) as usize];
                        d.entry(rep_len).or_insert_with(|| vec![]).push(rep_pos);
                    }
                }
            }
        }
    }

    let mut clustered = HashMap::new();
    for (rep_len, len_group) in d.into_iter() {
        let l = clustered.entry(rep_len).or_insert_with(|| vec![]);
        for cluster in connected_components(&to_graph::<(), ()>(&len_group)) {
            if cluster.len() > 1 {
                let mut cluster_l = cluster
                    .into_iter()
                    .map(|n| n.index())
                    .collect::<Vec<usize>>();
                cluster_l.sort();
                l.push(cluster_l);
            }
        }
    }
    Some(clustered)
}

#[pyfunction]
fn find_hairpins(sseq: String) -> PyResult<HashMap<usize, Vec<Vec<usize>>>> {
    match rust_find_hairpins(sseq) {
        None => Err(PyErr::new::<PyTypeError, _>("foo")),
        Some(a) => Ok(a),
    }
}

fn rust_find_hairpins(sseq: String) -> Option<HashMap<usize, Vec<Vec<usize>>>> {
    // move sseq to vector u8
    let seq: Vec<u8> = sseq.bytes().collect();
    let seq_rc = alphabets::dna::revcomp(&seq);

    // variables
    // let threshold_stem: usize = 8;
    // let threshold_stem: usize = 5;//2.25.25 settings trying to find smallest hairpins of stem len 2 - Drew
    // let threshold_stem: usize = 2;//2.25.25 settings trying to find smallest hairpins of stem len 2 - Drew
    let threshold_stem: usize = 4;//3.3.25 settings trying to find smallest hairpins of stem len 3 - Drew

    if seq.len() < threshold_stem * 2 {
        return None;
    }

    // create array
    let mut m = Matrix::from_elem((seq.len(), seq.len()), 0);

    // fill first row of array
    let first = seq_rc.iter().next().unwrap();
    for (i, b) in seq.iter().enumerate() {
        if b == first {
            m[[0, i]] = 1;
        }
    }

    // fill first column of array
    let first = seq.iter().next().unwrap();
    for (i, b) in seq_rc.iter().enumerate() {
        if b == first {
            m[[i, 0]] = 1;
        }
    }

    // fill above diagonal of the array, upper left
    for i in 1..seq.len() {
        for j in 1..seq.len() - i {
            if seq_rc[i] == seq[j] {
                m[[i, j]] = m[[i - 1, j - 1]] + 1
            }
        }
    }

    let mut d: HashMap<usize, Vec<Vec<usize>>> = HashMap::new();

    // search for hairpin patterns
    for i in 4..seq.len() - 4 {
        for j in 4..seq.len() - i {
            //if greater than stem length limit
            if (usize::from(m[[i, j]]) >= threshold_stem && m[[i, j]] > m[[i + 1, j + 1]])
                || (usize::from(m[[i, j]]) >= threshold_stem - 3
                    && m[[i, j]] > m[[i + 1, j + 1]]
                    && j + 1 == seq.len() - i - 1)
            {
                let a: usize = m[[i, j]].into();

                // if a local maximum (consolidates hairpins made of short repeats into the single largest hairpin)
                if m[[i, j]] > *m.slice(s![i, j - (a - 1)..j]).iter().max().unwrap()
                    && m[[i, j]] > *m.slice(s![i - (a - 1)..i, j]).iter().max().unwrap()
                {
                    let stem_len = m[[i, j]] as usize;
                    let stem_pos =
                        vec![(j - (stem_len - 1)) as usize, (seq.len() - i - 1) as usize];
                    d.entry(stem_len).or_insert_with(|| vec![]).push(stem_pos);
                }
            }
        }
    }
    Some(d)
}

#[pyfunction]
fn set_match_array(sseq: String, revcomp: bool, pm: &PyArray2<u8>) -> PyResult<()> {
    let seq: Vec<u8> = sseq.bytes().collect();
    if pm.shape() != vec![seq.len(), seq.len()] {
        return Err(PyErr::new::<PyTypeError, _>(format!(
            "expected numpy array w/ dim [{}, {}], got [{:?}]",
            seq.len(),
            seq.len(),
            pm.shape()
        )));
    }
    let mut m = pm.readwrite();

    if revcomp {
        let rc = alphabets::dna::revcomp(&seq);
        for i in 0..seq.len() {
            for j in 0..seq.len() {
                if seq[i] == rc[j] {
                    *m.get_mut([i, j]).unwrap() = 1;
                }
            }
        }
    } else {
        for i in 0..seq.len() {
            for j in 0..seq.len() {
                if seq[i] == seq[j] {
                    *m.get_mut([i, j]).unwrap() = 1;
                }
            }
        }
    }
    Ok(())
}

#[pyfunction]
fn gen_match_array(sseq: String, squery: String, pm: &PyArray2<u8>) -> PyResult<()> {
    let seq: Vec<u8> = sseq.bytes().collect();
    let query: Vec<u8> = squery.bytes().collect();
    if pm.shape() != vec![seq.len(), query.len()] {
        return Err(PyErr::new::<PyTypeError, _>(format!(
            "expected numpy array w/ dim [{}, {}], got [{:?}]",
            seq.len(),
            query.len(),
            pm.shape()
        )));
    }
    let mut m = pm.readwrite();

    for i in 0..seq.len() {
        for j in 0..query.len() {
            if seq[i] == query[j] {
                *m.get_mut([i, j]).unwrap() = 1;
            }
        }
    }
    Ok(())
}

// find near repeats
#[pyfunction]
fn find_near_repeats(sequence: &str) -> PyResult<Vec<HashMap<&str, i32>>> {
    let seed = 6;
    let gap_max = 2;
    // let min_length = 10;
    let min_length = 8;// 3.4.25
    let window = 10;
    let tolerance = 2;
    let repeat_log = near_repeat_log(
        sequence,
        seed,
        gap_max,
        mismatch_tolerance,
        min_length,
        window,
        tolerance,
    );
    match repeat_log_to_list(repeat_log) {
        None => Err(PyErr::new::<PyTypeError, _>("foo")),
        Some(m) => Ok(m),
    }
}

fn initialize_near_repeat_array(sequence: &str, features: usize) -> Vec<Vec<Vec<i32>>> {
    // returns len x len array
    // each element of returned array is a list
    // features of list include:
    // 0 lowest index of repeat start
    // 1 highest index of repeat start
    // 2 number of matches
    // 3 number of mismatches
    let n = sequence.len();
    let unit = vec![0; features];
    let row = vec![unit.clone(); n];
    vec![row.clone(); n]
}

fn mismatch_tolerance(match_count: i32, mismatch_count: i32, window: i32, tolerance: i32) -> bool {
    // determine if a position is within near repeat limits
    if (match_count / window) + 1 >= mismatch_count / tolerance {
        return true;
    } else {
        return false;
    }
}

fn near_repeat_log(
    sequence: &str,
    seed: i32,
    gap_max: usize,
    score_func: fn(i32, i32, i32, i32) -> bool,
    min_length: usize,
    window: i32,
    tolerance: i32,
) -> HashMap<i32, HashMap<i32, HashMap<&str, i32>>> {
    // different types of sequence representation
    let ssequence: String = sequence.to_string();
    let isequence: Vec<u8> = ssequence.bytes().collect();

    // data structures for finding and reporting repeats
    let mut R = initialize_near_repeat_array(sequence, 4);
    let mut repeat_log: HashMap<i32, HashMap<i32, HashMap<&str, i32>>> = HashMap::new();

    // first row of array
    for j in 1..sequence.len() {
        if isequence[0] == isequence[j] {
            R[0][j][1] = j as i32;
            R[0][j][2] = 1;
        }
    }

    // fill in top right triangle of array
    for i in 1..sequence.len() {
        for j in i + min_length..sequence.len() {
            // create i32 types for indexes i and j
            let I: i32 = i.try_into().unwrap();
            let J: i32 = j.try_into().unwrap();
            if isequence[i] == isequence[j] {
                // find max match count within allowed gap distance
                let field_dim = cmp::min(gap_max + 1, i);
                // if position is far enough from an edge to consider mismatches
                if field_dim > 0 {
                    let candidates: Vec<_> = R[i - field_dim..i]
                        .iter()
                        .flat_map(|row| row[j - field_dim..j].iter())
                        .collect();
                    // sort and select top candidate
                    let best_repeat = {
                        let mut sorted_candidates = candidates.clone();
                        sorted_candidates.sort_by(|a, b| {
                            let a_tuple = (a[2], -a[3], a[0], a[1]);
                            let b_tuple = (b[2], -b[3], b[0], b[1]);
                            a_tuple.cmp(&b_tuple)
                        });
                        sorted_candidates.reverse();
                        sorted_candidates[0].clone()
                    };
                    // allow mismatches only if repeat is longer than seed
                    if best_repeat[2] >= seed {
                        for m in 0..field_dim {
                            for n in 0..field_dim {
                                if R[i - field_dim + m][j - field_dim + n] == best_repeat {
                                    // update repeat start indexes
                                    R[i][j][0] = best_repeat[0];
                                    R[i][j][1] = best_repeat[1];
                                    // update match base count
                                    R[i][j][2] = best_repeat[2] + 1;
                                    // update mismatch count
                                    R[i][j][3] = best_repeat[3]
                                        + cmp::max(field_dim - m - 1, field_dim - n - 1) as i32;
                                    if score_func(R[i][j][2], R[i][j][3], window, tolerance) == true
                                    {
                                        // log if long enough, not overlapping, max length of local repeats
                                        if I - R[i][j][0] + 1 >= min_length.try_into().unwrap()
                                            && R[i][j][1] - R[i][j][0] >= I - R[i][j][0] + 1
                                            && R[i][j][2]
                                                >= (1..R[i][j][2] as usize)
                                                    .map(|x| R[i][j - x][2])
                                                    .max()
                                                    .unwrap()
                                            && R[i][j][2]
                                                >= (1..R[i][j][2] as usize)
                                                    .map(|x| R[i - x][j][2])
                                                    .max()
                                                    .unwrap()
                                        {
                                            let start1 = R[i][j][0];
                                            let start2 = R[i][j][1];
                                            if !repeat_log.contains_key(&start1) {
                                                repeat_log.insert(
                                                    start1,
                                                    std::collections::HashMap::new(),
                                                );
                                            }

                                            if !repeat_log[&start1].contains_key(&start2) {
                                                repeat_log.get_mut(&start1).unwrap().insert(
                                                    start2,
                                                    std::collections::HashMap::from([
                                                        ("length", 0),
                                                        ("match_count", 0),
                                                        ("mismatch_count", 0),
                                                    ]),
                                                );
                                            }
                                            let values = repeat_log
                                                .get_mut(&start1)
                                                .unwrap()
                                                .get_mut(&start2)
                                                .unwrap();
                                            values.insert("length", I - start1 + 1);
                                            values.insert("match_count", R[i][j][2]);
                                            values.insert("mismatch_count", R[i][j][3]);
                                        }
                                    }
                                    // if mismatch too high, reset element to zeros and break repeat
                                    else {
                                        R[i][j] = R[0][0].clone();
                                    } // assumes first element of array remains zero throughout
                                }
                            }
                        }
                    }
                    // if repeat length is less than seed length, extend along diagonal
                    else {
                        R[i][j][2] = R[i - 1][j - 1][2] + 1;
                        R[i][j][0] = I - R[i][j][2] + 1;
                        R[i][j][1] = J - R[i][j][2] + 1;
                    }
                } else {
                    R[i][j][2] = R[i - 1][j - 1][2] + 1;
                    R[i][j][0] = I - R[i][j][2] + 1;
                    R[i][j][1] = J - R[i][j][2] + 1;
                }
            }
        }
    }
    repeat_log
}

fn repeat_log_to_list(
    repeat_log: HashMap<i32, HashMap<i32, HashMap<&str, i32>>>,
) -> Option<Vec<HashMap<&str, i32>>> {
    // collapse to longest repeat
    // find sets of repeats that overlap each other and select the longest to represent
    let mut repeat_list: Vec<HashMap<&str, i32>> = Vec::new();
    for (start1, value1) in &repeat_log {
        let mut start2_list: Vec<Vec<i32>> = Vec::new();
        for start2 in value1.keys() {
            let shared: Vec<i32> = value1
                .keys()
                .filter(|&&x| (start2 - x).abs() <= value1[&start2]["length"])
                .cloned()
                .collect();
            if start2_list.is_empty() {
                start2_list.push(shared);
            } else {
                let mut is_shared = Cell::new(false);
                for x in &mut start2_list {
                    if shared.iter().any(|v| x.contains(v)) {
                        x.extend(shared.clone());
                        is_shared = true.into();
                    }
                }
                if !is_shared.get() {
                    start2_list.push(shared);
                }
            }
        }
        let mut start2_sets: Vec<Vec<i32>> = Vec::new();
        for mut x in start2_list {
            x.sort_unstable();
            x.dedup();
            start2_sets.push(x);
        }
        for s in start2_sets {
            let longest = s
                .iter()
                .map(|x| {
                    let length = repeat_log[&start1][*&x]["length"];
                    let match_count = repeat_log[&start1][*&x]["match_count"];
                    let mismatch_count = repeat_log[&start1][*&x]["mismatch_count"];
                    (length, *x, match_count, mismatch_count)
                })
                .max_by_key(|tuple| *tuple)
                .unwrap();
            repeat_list.push(HashMap::from([
                ("position1", start1.clone()),
                ("position2", longest.1),
                ("length", longest.0),
                ("match_count", longest.2),
                ("mismatch_count", longest.3),
            ]));
        }
    }
    Some(repeat_list)
}

// find near hairpins
#[pyfunction]
fn find_near_hairpins(sequence: &str) -> PyResult<Vec<HashMap<&str, i32>>> {
    let seed = 6;
    let gap_max = 2;
    // let min_length = 10;
    // let min_length = 5;//2.25.25 settings to find all hairpins
    let min_length = 8;//3.3.25 settings to find all hairpins
    let window = 10;
    let tolerance = 2;
    let hp_log = near_hairpin_log(
        sequence,
        seed,
        gap_max,
        mismatch_tolerance,
        min_length,
        window,
        tolerance,
    );
    match hp_log_to_list(sequence, hp_log) {
        None => Err(PyErr::new::<PyTypeError, _>("foo")),
        Some(m) => Ok(m),
    }
}
fn near_hairpin_log(
    sequence: &str,
    seed: i32,
    gap_max: usize,
    score_func: fn(i32, i32, i32, i32) -> bool,
    min_length: usize,
    window: i32,
    tolerance: i32,
) -> HashMap<i32, HashMap<i32, HashMap<&str, i32>>> {
    // different sequence and reverse complement representations
    
    let ssequence: String = sequence.to_string();
    let isequence: Vec<u8> = ssequence.bytes().collect();
	
	let isequence_rc = alphabets::dna::revcomp(isequence.clone());

    // data structures for finding and reporting repeats
    let mut R = initialize_near_repeat_array(sequence, 4);
    let mut hp_log: HashMap<i32, HashMap<i32, HashMap<&str, i32>>> = HashMap::new();

    // first row of array
    for j in 0..sequence.len() {
        if isequence[0] == isequence_rc[j] {
            R[0][j][1] = j as i32;
            R[0][j][2] = 1;
        }
    }
	// first column of array
    for i in 0..sequence.len() {
        if isequence[i] == isequence_rc[0] {
            R[i][0][1] = i as i32;
            R[i][0][2] = 1;
        }
	}
	
    // fill in top left triangle of array
    for i in 1..sequence.len() {
        for j in 1..sequence.len() - i {
            // create i32 types for indexes i and j
            let I: i32 = i.try_into().unwrap();
            let J: i32 = j.try_into().unwrap();
			//println!("(I): {:#?}", I);
			//println!("(J): {:#?}", J);
            if isequence[i] == isequence_rc[j] {
				//println!("i: {i} j: {j} i char: {} j char: {}", isequence[i] as char, isequence_rc[j] as char);
                // find max match count within allowed gap distance
                let field_dim = cmp::min(gap_max + 1, cmp::min(i, j));
                // if position is far enough from an edge to consider mismatches
                if field_dim > 0 {
                    let candidates: Vec<_> = R[i - field_dim..i]
                        .iter()
                        .flat_map(|row| row[j - field_dim..j].iter())
                        .collect();
                    // sort and select top candidate
                    let best_repeat = {
                        let mut sorted_candidates = candidates.clone();
                        sorted_candidates.sort_by(|a, b| {
                            let a_tuple = (a[2], -a[3], a[0], a[1]);
                            let b_tuple = (b[2], -b[3], b[0], b[1]);
                            a_tuple.cmp(&b_tuple)
                        });
                        sorted_candidates.reverse();
                        sorted_candidates[0].clone()
                    };
                    // allow mismatches only if repeat is longer than seed
                    if best_repeat[2] >= seed {
                        for m in 0..field_dim {
                            for n in 0..field_dim {
                                if R[i - field_dim + m][j - field_dim + n] == best_repeat {
                                    // update repeat start indexes
                                    R[i][j][0] = best_repeat[0];
                                    R[i][j][1] = best_repeat[1];
                                    // update match base count
                                    R[i][j][2] = best_repeat[2] + 1;
                                    // update mismatch count
                                    R[i][j][3] = best_repeat[3]
                                        + cmp::max(field_dim - m - 1, field_dim - n - 1) as i32;
                                    if score_func(R[i][j][2], R[i][j][3], window, tolerance) {
                                        // log if long enough, not overlapping, max length of local repeats
                                        if I - R[i][j][0] + 1 >= min_length.try_into().unwrap()
                                            && R[i][j][2]
                                                > (1..R[i][j][2] as usize)
                                                    .map(|x| R[i][j - x][2])
                                                    .max()
                                                    .unwrap()
                                            && R[i][j][2]
                                                > (1..R[i][j][2] as usize)
                                                    .map(|x| R[i - x][j][2])
                                                    .max()
                                                    .unwrap()
                                        {
                                            let start1 = R[i][j][0];
                                            let start2 = R[i][j][1];
                                            if !hp_log.contains_key(&start1) {
                                                hp_log.insert(
                                                    start1,
                                                    std::collections::HashMap::new(),
                                                );
                                            }

                                            if !hp_log[&start1].contains_key(&start2) {
                                                hp_log.get_mut(&start1).unwrap().insert(
                                                    start2,
                                                    std::collections::HashMap::from([
                                                        ("length", 0),
                                                        ("match_count", 0),
                                                        ("mismatch_count", 0),
                                                    ]),
                                                );
                                            }
                                            let values = hp_log
                                                .get_mut(&start1)
                                                .unwrap()
                                                .get_mut(&start2)
                                                .unwrap();
                                            values.insert("length", I - start1 + 1);
                                            values.insert("match_count", R[i][j][2]);
                                            values.insert("mismatch_count", R[i][j][3]);
                                        }
                                    }
                                    // if mismatch too high, reset element to zeros and break repeat
                                    else {
                                        R[i][j] = R[sequence.len()-1][sequence.len()-1].clone();
                                    } // assumes bottom corner element of array remains zero throughout
                                }
                            }
                        }
                    }
                    // if repeat length is less than seed length, extend along diagonal
                    else {
                        R[i][j][2] = R[i - 1][j - 1][2] + 1;
                        R[i][j][0] = I - R[i][j][2] + 1;
                        R[i][j][1] = J - R[i][j][2] + 1;
                    }
                }
                // if position prevents considering mismatches
                else {
                    R[i][j][2] = R[i - 1][j - 1][2] + 1;
                    R[i][j][0] = I - R[i][j][2] + 1;
                    R[i][j][1] = J - R[i][j][2] + 1;
                }
            }
        }
    }
    hp_log
}

fn hp_log_to_list<'a>(
    sequence: &'a str,
    hp_log: HashMap<i32, HashMap<i32, HashMap<&'a str, i32>>>,
) -> Option<Vec<HashMap<&'a str, i32>>> {
    // collapse to longest repeat
    // find sets of repeats that overlap each other and select the longest to represent
    let mut hp_list: Vec<HashMap<&str, i32>> = Vec::new();
    for (start1, value1) in &hp_log {
        let mut start2_list: Vec<Vec<i32>> = Vec::new();
        for start2 in value1.keys() {
            let shared: Vec<i32> = value1
                .keys()
                .filter(|&&x| (start2 - x).abs() <= value1[&start2]["length"])
                .cloned()
                .collect();
            if start2_list.is_empty() {
                start2_list.push(shared);
            } else {
                let mut is_shared = Cell::new(false);
                for x in &mut start2_list {
                    if shared.iter().any(|v| x.contains(v)) {
                        x.extend(shared.clone());
                        is_shared = true.into();
                    }
                }
                if !is_shared.get() {
                    start2_list.push(shared);
                }
            }
        }
        let mut start2_sets: Vec<Vec<i32>> = Vec::new();
        for mut x in start2_list {
            x.sort_unstable();
            x.dedup();
            start2_sets.push(x);
        }
        for s in start2_sets {
            let longest = s
                .iter()
                .map(|x| {
                    let length = hp_log[&start1][*&x]["length"];
                    let match_count = hp_log[&start1][*&x]["match_count"];
                    let mismatch_count = hp_log[&start1][*&x]["mismatch_count"];
                    (length, *x, match_count, mismatch_count)
                })
                .max_by_key(|tuple| *tuple)
                .unwrap();
            hp_list.push(HashMap::from([
                ("position1", start1.clone()),
                ("position2", sequence.len() as i32 - longest.1 - longest.0),
                ("length", longest.0),
                ("match_count", longest.2),
                ("mismatch_count", longest.3),
            ]));
        }
    }
    Some(hp_list)
}

type HS = [f64; 2];
const D_H: usize = 0;
const D_S: usize = 1;

fn hs_lookup(first: u8, second: u8) -> Option<HS> {
    match (first, second) {
        (b'A', b'A') | (b'T', b'T') => Some([-7.9, -22.2]),
        (b'A', b'T') => Some([-7.2, -20.4]),
        (b'T', b'A') => Some([-7.2, -21.3]),
        (b'C', b'A') | (b'T', b'G') => Some([-8.5, -22.7]),
        (b'G', b'T') | (b'A', b'C') => Some([-8.4, -22.4]),
        (b'C', b'T') | (b'A', b'G') => Some([-7.8, -21.0]),
        (b'G', b'A') | (b'T', b'C') => Some([-8.2, -22.2]),
        (b'C', b'G') => Some([-10.6, -27.2]),
        (b'G', b'C') => Some([-9.8, -24.4]),
        (b'G', b'G') | (b'C', b'C') => Some([-8.0, -19.9]),
        _ => None,
    }
}

#[pyfunction]
fn tm_nn(sseq: String) -> f64 {
    /*
    minimal port of Biopython's Tm_NN
    doesn't cover: dangling ends, mismatches (uses reverse-comp)
    DNA_NN3 = {
        "init": (0, 0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8),
        "init_oneG/C": (0, 0), "init_allA/T": (0, 0), "init_5T/A": (0, 0),
        "sym": (0, -1.4),
        "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
        "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
        "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9)}
         */

    let seq: Vec<u8> = sseq.to_uppercase().bytes().collect();

    let nn3_init: HS = [0.0, 0.0];
    let nn3_initAT: HS = [2.3, 4.1];
    let nn3_initGC: HS = [0.1, -2.8];
    let _nn3_sym: HS = [0.0, -1.4];
    let dnac1: f64 = 500.0;
    let dnac2: f64 = 0.0;
    let mut delta_H = nn3_init[D_H];
    let mut delta_S = nn3_init[D_S];

    let end_AT = if seq[0] == b'A' || seq[0] == b'T' {
        1.0
    } else {
        0.0
    } + if seq[seq.len() - 1] == b'A' || seq[seq.len() - 1] == b'T' {
        1.0
    } else {
        0.0
    };
    delta_H += end_AT * nn3_initAT[D_H];
    delta_S += end_AT * nn3_initAT[D_S];
    delta_H += (2.0 - end_AT) * nn3_initGC[D_H];
    delta_S += (2.0 - end_AT) * nn3_initGC[D_S];
    for p in seq.windows(2) {
        let lk =
            hs_lookup(p[0], p[1]).expect(&format!("unknown: {}{}", p[0] as char, p[1] as char));
        delta_H += lk[D_H];
        delta_S += lk[D_S];
    }
    let k = (dnac1 - (dnac2 / 2.0)) * 1e-9;
    let R = 1.987; // universal gas constant in Cal/degrees C*Mol

    let corr = {
        let na: f64 = 100.0;
        let mg: f64 = 3.0;
        let k: f64 = 0.0;
        let tris: f64 = 0.0;
        let _dNTPs: f64 = 0.0;

        let mut mon_: f64 = na + k + tris / 2.0;
        let _mg_: f64 = mg * 1e-3;
        mon_ += 120.0 * mg.sqrt();
        let mon: f64 = mon_ * 1e-3;
        16.6 * mon.log10()
    };
    corr + (1000.0 * delta_H) / (delta_S + (R * k.ln())) - 273.15
}

fn compute_hashes(seq: &[u8], k: usize) -> Vec<usize> {
    let md: usize = 1e9 as usize + 7;
    let mut most_significant: usize = 1;
    for _i in 1..k {
        most_significant = most_significant * 4 % md;
    }
    fn nucval(base: u8) -> usize {
        match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        }
    }

    let extendr = |prev, nuc| ((prev * 4) + nucval(nuc)) % md;
    let remove_left = |prev, nuc| (prev + md - (most_significant * nucval(nuc)) % md) % md;
    let _h: usize = 0;
    let h = seq[..k].iter().fold(0, |acc, nuc| extendr(acc, *nuc));
    let mut hashes = vec![h];
    let mut curr = h;
    hashes.extend((1..seq.len() - k + 1).map(|pos| {
        curr = remove_left(curr, seq[pos - 1]);
        curr = extendr(curr, seq[pos + k - 1]);
        curr
    }));
    hashes.extend((seq.len() - k + 1..seq.len()).map(|_| 0));

    hashes
}

#[pyfunction]
fn find_repeats_to_separate(sseq: String) -> PyResult<HashMap<usize, Vec<(usize, usize)>>> {
    match rust_find_repeats_to_separate(sseq) {
        None => Err(PyErr::new::<PyTypeError, _>("foo")),
        Some(m) => Ok(m),
    }
}

fn get_closest_repeat(seq: &[u8], h: &[usize], i: usize, k: usize) -> Option<(usize, usize)> {
    for j in i + k..seq.len() - k + 1 {
        if h[i] == h[j] && seq[i..i + k] == seq[j..j + k] {
            let mut l = k;
            while j + l < seq.len() && seq[i + l] == seq[j + l] {
                l += 1;
            }
            let hs: HashSet<u8> = seq[i..j + l].iter().cloned().collect();
            if i + l <= j && hs.len() > 1 {
                // no overlap and not within the same homopolymer
                return Some((j, l));
            }
        }
    }
    return None;
}

fn get_closest_rc_repeat(
    seq: &[u8],
    h: &[usize],
    h_rc: &[usize],
    i: usize,
    k: usize,
) -> Option<(usize, usize)> {
    for j in i + 2 * k..seq.len() {
        if h[i] == h_rc[j]
            && seq[i..i + k] == alphabets::dna::revcomp(seq[j - k + 1..j + 1].to_vec())
        {
            let mut l = k;
            while i + l < j - l && seq[i + l] == alphabets::dna::complement(seq[j - l]) {
                l += 1;
            }
            return Some((j - l + 1, l));
        }
    }
    return None;
}

fn rust_find_repeats_to_separate(sseq: String) -> Option<HashMap<usize, Vec<(usize, usize)>>> {
    let seq: Vec<u8> = sseq.bytes().collect();
    let k: usize = 41; // repeat lower limit
    let mut d: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();

    if seq.len() < 2 * k {
        return Some(d);
    }

    let h = compute_hashes(&seq, k);
    let seq_rc: Vec<u8> = alphabets::dna::revcomp(&seq);
    let h_rc = compute_hashes(&seq_rc, k)
        .into_iter()
        .rev()
        .collect::<Vec<_>>();

    let mut last_match: usize = 0;
    for i in 0..seq.len() - (2 * k) {
        if let Some((j, l)) = get_closest_repeat(&seq, &h, i, k) {
            if last_match != j - 1 {
                d.entry(l).or_insert_with(|| vec![]).push((i, j));
            }
            last_match = j;
        } else {
            // no repeat found
            last_match = 0;
        }
    }

    last_match = 0;
    for i in 0..seq.len() - (2 * k) {
        if let Some((j, l)) = get_closest_rc_repeat(&seq, &h, &h_rc, i, k) {
            if last_match != j {
                d.entry(l).or_insert_with(|| vec![]).push((i, j));
            }
            last_match = j;
        } else {
            // no revcomp repeat found
            last_match = 0;
        }
    }

    return Some(d);
}

fn ops_to_cigarstring<'a>(ops: std::slice::Iter<'a, Op>) -> String {
    let mut cigar = String::new();
    let mut count: usize = 0;
    let mut prev = Op::Match;
    for (idx, op) in ops.enumerate() {
        if op == &prev {
            count += 1;
        } else {
            if idx > 0 {
                cigar.push_str(&count.to_string());
                cigar.push(match prev {
                    Op::Match => '=',
                    Op::Subst => 'X',
                    Op::Ins => 'I',
                    Op::Del => 'D',
                    _ => unreachable!(),
                });
            }
            prev = op.clone();
            count = 1;
        }
    }
    cigar.push_str(&count.to_string());
    cigar.push(match prev {
        Op::Match => '=',
        Op::Subst => 'X',
        Op::Ins => 'I',
        Op::Del => 'D',
        _ => unreachable!(),
    });
    cigar
    //cigar.chars().rev().collect()
}

fn rust_pairwise_local(
    subj: &[u8],
    query: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> (i32, usize, usize, String) {
    let mut aligner = pairwise::Aligner::with_capacity(
        subj.len(),
        query.len(),
        gap_open,
        gap_extend,
        |a: u8, b: u8| if a == b { match_score } else { mismatch_score },
    );
    let aln = aligner.local(&query, &subj);
    (
        aln.score,
        aln.xstart,
        aln.ystart,
        // aln.cigar(false), // can't use!
        ops_to_cigarstring(aln.operations.iter()),
    )
}

#[pyfunction]
fn pairwise_local(
    subj: String,
    query: String,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<(i32, usize, usize, String)> {
    Ok(rust_pairwise_local(
        subj.as_bytes(),
        query.as_bytes(),
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    ))
}

/// A Python module implemented in Rust.
#[pymodule]
fn speedy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(find_repeats, m)?)?;
    m.add_function(wrap_pyfunction!(find_near_repeats, m)?)?;
    m.add_function(wrap_pyfunction!(find_hairpins, m)?)?;
    m.add_function(wrap_pyfunction!(find_near_hairpins, m)?)?;
    m.add_function(wrap_pyfunction!(set_match_array, m)?)?;
    m.add_function(wrap_pyfunction!(gen_match_array, m)?)?;
    m.add_function(wrap_pyfunction!(tm_nn, m)?)?;
    m.add_function(wrap_pyfunction!(find_repeats_to_separate, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_local, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_str2() {
        assert_eq!(
            rust_find_repeats(String::from("ATATATATATATATATATATATAT"))
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn test_rep_to_separate() {
        let rep = "ATCCTTGAAGGAGGAATTATTAGCAGGTGCGATAATCGTATGG";
        let seq = format!("{}{}", rep, rep);
        let repeats = rust_find_repeats_to_separate(seq).unwrap();
        assert_eq!(repeats.len(), 1);
        assert_eq!(repeats[&43], [(0, 43)]);
    }

    #[test]
    fn test_rep_to_separate_revcomp() {
        let rep = "ATCCTTGAAGGAGGAATTATTAGCAGGTGCGATAATCGTATGG";
        let rc = String::from_utf8(alphabets::dna::revcomp(rep.as_bytes())).unwrap();
        let seq = format!("{}{}", rep, rc);
        let repeats = rust_find_repeats_to_separate(seq).unwrap();
        assert_eq!(repeats.len(), 1);
        println!("{:?}", repeats);
        assert_eq!(repeats[&43], [(0, 43)]);
    }

    #[test]
    fn test_rep_to_separate_homopolymers() {
        let hp1 = "A".repeat(50);
        let hp2 = "C".repeat(100);
        let seq = format!("{}{}{}", hp1, hp2, hp1);
        let repeats = rust_find_repeats_to_separate(seq).unwrap();
        assert_eq!(repeats.len(), 10);
        assert_eq!(repeats[&50], [(0, 150)]);
    }

    #[test]
    fn test_rep_to_separate_overlapping_repeat() {
        let rep1 = "CACACCTTATAT";
        let rep2 = "CGAACTTAAGCCATCCGATAAGACAGTCC";
        let seq = format!("{}{}{}{}{}{}{}", rep1, rep2, rep1, rep2, rep1, rep2, rep1);
        let repeats = rust_find_repeats_to_separate(seq).unwrap();
        assert_eq!(repeats.len(), 1);
        assert_eq!(repeats[&53], [(0, 82)]);
    }

    #[test]
    fn test_rep_to_separate_subrepeat() {
        let subrep = "ACCATTCAACACATGGAACGCTCTTTATCACCTCCCAACAT";
        let rep = format!("GGTC{}", subrep);
        let seq = format!("{}TT{}AA{}", rep, subrep, rep);
        let repeats = rust_find_repeats_to_separate(seq).unwrap();
        assert_eq!(repeats.len(), 2);
        assert_eq!(repeats[&45], [(0, 90)]);
        assert_eq!(repeats[&41], [(4, 47), (47, 94)]);
    }

    #[test]
    fn test_pw_simple() {
        assert_eq!(
            rust_pairwise_local("ACGT".as_bytes(), "ACGT".as_bytes(), 1, -1, -1, -1),
            (4, 0, 0, "4=".to_string())
        );
    }

    #[test]
    fn test_pw_ins() {
        assert_eq!(
            rust_pairwise_local("ATGCATGC".as_bytes(), "ATGCTATGC".as_bytes(), 1, -1, -1, -1),
            (6, 0, 0, "4=1I4=".to_string())
        );
    }

    #[test]
    fn test_pw_del() {
        assert_eq!(
            rust_pairwise_local("ATGCATGC".as_bytes(), "ATGATGC".as_bytes(), 1, -1, -1, -1),
            (5, 0, 0, "3=1D4=".to_string())
        );
    }

    #[test]
    fn test_pw_last() {
        assert_eq!(
            rust_pairwise_local("ATGCATGC".as_bytes(), "ATGCATC".as_bytes(), 2, -1, -1, -1),
            (12, 0, 0, "6=1D1=".to_string())
        );
    }
}
