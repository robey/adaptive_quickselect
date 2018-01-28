#[cfg(test)]
extern crate rand;

use std::fmt::Debug;

/// Returns the index of the median of r[a], r[b], r[c].
#[inline]
fn median_index<A: PartialOrd>(r: &[A], a: usize, b: usize, c: usize) -> usize {
  if r[a] > r[c] {
    if r[b] > r[a] { a } else { if r[b] < r[c] { c } else { b } }
  } else {
    if r[b] > r[c] { c } else { if r[b] < r[a] { a } else { b } }
  }
}

/// Tukey's ninther: compute the median of each triplet at a time, then the
/// median of those medians, returning that index.
#[inline]
fn ninther<A: PartialOrd>(
  r: &mut [A],
  a: usize, b: usize, c: usize,
  d: usize, e: usize, f: usize,
  g: usize, h: usize, i: usize
) -> usize {
  median_index(r, median_index(r, a, b, c), median_index(r, d, e, f), median_index(r, g, h, i))
}

/// Given a mostly-partitioned array with an unpartitioned gap on the right
/// side of the pivot in [hi, right), move those into place, shifting the
/// pivot right as necessary. Return the new pivot index.
/// [ p > > > > x x x x > > ]  ->  [ < < < < p > > > > > > ]
///   pivot     hi      right                pivot     right
///
fn expand_partition_right<A: PartialOrd>(r: &mut [A], pivot: usize, hi: usize, right: usize) -> usize {
  let mut p = pivot;
  for i in hi..right {
    if r[i] < r[pivot] {
      p += 1;
      r.swap(p, i);
    }
  }
  r.swap(p, pivot);
  p
}

/// Same as `expand_partition_right`, but in reverse.
/// [ < < x x x x < < < p ]
///       left    lo    pivot
fn expand_partition_left<A: PartialOrd>(r: &mut [A], pivot: usize, lo: usize, left: usize) -> usize {
  let mut p = pivot;
  for i in (left..lo).rev() {
    if r[i] > r[pivot] {
      p -= 1;
      r.swap(p, i);
    }
  }
  r.swap(p, pivot);
  p
}

/// Given a pivot and a [lo, hi) span around it which is already partitioned,
/// partition the rest, returning the new pivot position.
fn expand_partition<A: PartialOrd + Debug>(r: &mut [A], lo: usize, pivot: usize, hi: usize) -> usize {
  debug_assert!(lo <= pivot && pivot < hi && hi <= r.len());
  let mut left = 0;
  let mut right = r.len() - 1;

  loop {
    while left < lo && r[left] <= r[pivot] { left += 1 }
    if left == lo { return expand_partition_right(r, pivot, hi, right + 1); }
    while right >= hi && r[right] >= r[pivot] { right -= 1 }
    if right < hi { return expand_partition_left(r, pivot, lo, left); }

    r.swap(left, right);
    left += 1;
    right -= 1;
  }
}

/// partition by the element at `k`, and return the new position of that element.
fn partition_hoare<A: PartialOrd>(r: &mut [A], k: usize) -> usize {
  debug_assert!(k < r.len());
  r.swap(0, k);
  let mut lo = 1;
  let mut hi = r.len() - 1;

  while lo <= hi {
    while lo <= hi && r[lo] < r[0] { lo += 1 }
    while lo <= hi && r[hi] > r[0] { hi -= 1 }
    if lo <= hi {
      r.swap(lo, hi);
      lo += 1;
      hi -= 1;
    }
  }

  lo -= 1;
  r.swap(lo, 0);
  lo
}

/// Partition by finding a pivot that approximates the median. This algorithm
/// is described in great detail in the paper (listed in the README). It
/// samples 3 subsets of the array, from 3 equidistant ranges, and computes
/// the median-of-median (`ninther`) of 9 items at a time, in groups of 3.
/// These medians are moved into the center of the array, partitioned around
/// the center, and then the partioning is expanded to cover the whole array.
fn partition_ninthers<A: PartialOrd + Debug>(r: &mut [A]) -> usize {
  debug_assert!(r.len() >= 12);
  let frac = if r.len() <= 1024 {
    r.len() / 12
  } else {
    if r.len() <= 128 * 1024 { r.len() / 64 } else { r.len() / 1024 }
  };

  let pivot = frac / 2;
  let lo = r.len() / 2 - pivot;
  let hi = lo + frac;

  debug_assert!(lo >= frac * 4);
  debug_assert!(r.len() - hi >= frac * 4);
  debug_assert!(lo / 2 >= pivot);

  let gap = (r.len() - 9 * frac) / 4;
  let mut a = lo - 4 * frac - gap;
  let mut b = hi + gap;
  for i in lo..hi {
    let k = ninther(r, a, i - frac, b, a + 1, i, b + 1, a + 2, i + frac, b + 2);
    r.swap(k, i);
    a += 3;
    b += 3;
  }

  adaptive_quickselect(&mut r[lo..hi], pivot);
  expand_partition(r, lo, lo + pivot, hi)
}

/// Walk the entire array in chunks, collecting the minimum of each chunk
/// into the first 2*k elements, then recursively find the median of those.
fn partition_minima<A: PartialOrd + Debug>(r: &mut [A], k: usize) -> usize {
  debug_assert!(r.len() >= 2);
  debug_assert!(k > 0 && k * 4 <= r.len());

  let subset = k * 2;
  let span = (r.len() - subset) / subset;
  let start = 0;
  let end = start + subset;
  debug_assert!(span > 0);

  let mut chunk = subset;
  for i in start..end {
    let mut index = chunk;
    for j in (chunk + 1)..(chunk + span) {
      if r[j] < r[index] { index = j }
    }
    if r[index] < r[i] { r.swap(index, i) }
    chunk += span;
  }
  debug_assert!(chunk <= r.len() && chunk + subset > r.len());

  adaptive_quickselect(&mut r[start..end], k);
  expand_partition(r, start, k, end)
}

/// Walk the entire array in chunks, collecting the maximum of each chunk
/// into the last 2*k elements, then recursively find the median of those.
fn partition_maxima<A: PartialOrd + Debug>(r: &mut [A], k: usize) -> usize {
  debug_assert!(r.len() >= 2);
  debug_assert!(k < r.len() && k * 4 >= r.len() * 3);

  let subset = (r.len() - k) * 2;
  let span = (r.len() - subset) / subset;
  let start = r.len() - subset;
  let end = r.len();
  debug_assert!(span > 0);

  let mut chunk = start - subset * span;
  for i in start..end {
    let mut index = chunk;
    for j in (chunk + 1)..(chunk + span) {
      if r[j] > r[index] { index = j }
    }
    if r[index] > r[i] { r.swap(index, i) }
    chunk += span;
  }
  println!("{}, {}, len={} k={}", chunk, subset, r.len(), k);
  debug_assert!(chunk == r.len() - subset);

  let pivot = r.len() - k;
  adaptive_quickselect(&mut r[start..end], pivot);
  expand_partition(r, start, k, end)
}

/// Select the `k`th item from a hypothetically-sorted version of `r`.
/// For example, the median of `r` is `adaptive_quickselect(r, r.len() / 2)`.
/// `r` will be modified by partially sorting (partitioning) around various
/// pivots until the `k`th item is in the right place. No guarantee is made
/// about any other position.
pub fn adaptive_quickselect<A: PartialOrd + Debug>(r: &mut [A], k: usize) {
  debug_assert!(k < r.len());
  let last = r.len() - 1;

  if k == 0 {
    // minimum.
    let mut pivot = 0;
    for i in 1..r.len() {
      if r[i] < r[pivot] { pivot = i }
    }
    r.swap(0, pivot);
    return;
  }

  if k == last {
    // maximum.
    let mut pivot = 0;
    for i in 1..r.len() {
      if r[i] > r[pivot] { pivot = i }
    }
    r.swap(last, pivot);
    return;
  }

  let pivot = if r.len() <= 16 {
    partition_hoare(r, k)
  } else if k * 6 <= r.len() {
    partition_minima(r, k)
  } else if k * 6 >= r.len() * 5 {
    partition_maxima(r, k)
  } else {
    partition_ninthers(r)
  };

  // FIXME how to mark this is a tail call to make sure the compiler gets it?
  if pivot == k { return }
  if pivot > k {
    adaptive_quickselect(&mut r[..pivot], k);
  } else {
    let start = pivot + 1;
    adaptive_quickselect(&mut r[start..], k - start);
  }
}


// ----- tests

#[cfg(test)]
mod tests {
  use adaptive_quickselect;
  use rand::{Rng, SeedableRng, StdRng};

  fn make_random_sequence(generator: &mut StdRng, n: usize) -> Vec<u32> {
    (0..n).map(|_| generator.next_u32()).collect::<Vec<u32>>()
  }

  #[test]
  fn test_min() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..100 {
      let mut r = make_random_sequence(&mut generator, 100);
      adaptive_quickselect(&mut r, 0);
      for n in r[1..].iter() {
        assert!(*n >= r[0]);
      }
    }
  }

  #[test]
  fn test_max() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..100 {
      let mut r = make_random_sequence(&mut generator, 100);
      adaptive_quickselect(&mut r, 99);
      for n in r[0..99].iter() {
        assert!(*n <= r[99]);
      }
    }
  }

  #[test]
  fn test_hoare() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..50000 {
      let mut r = make_random_sequence(&mut generator, 16);
      adaptive_quickselect(&mut r, 8);
      let answer = r[8];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[8], answer);
    }
  }

  #[test]
  fn test_ninthers_median_1k() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..5000 {
      let mut r = make_random_sequence(&mut generator, 1024);
      adaptive_quickselect(&mut r, 512);
      let answer = r[512];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[512], answer);
    }
  }

  #[test]
  fn test_ninthers_median_10k() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..100 {
      let mut r = make_random_sequence(&mut generator, 10240);
      adaptive_quickselect(&mut r, 5120);
      let answer = r[5120];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[5120], answer);
    }
  }

  #[test]
  fn test_ninthers_extrema_10k() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..100 {
      let mut r = make_random_sequence(&mut generator, 10240);
      adaptive_quickselect(&mut r, 50);
      let answer = r[50];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[50], answer);
    }

    for _ in 0..100 {
      let mut r = make_random_sequence(&mut generator, 10240);
      adaptive_quickselect(&mut r, 10201);
      let answer = r[10201];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[10201], answer);
    }
  }

  #[test]
  fn test_ninthers_any_1k() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..5000 {
      let mut r = make_random_sequence(&mut generator, 1024);
      let k = (generator.next_u32() & 0x3ff) as usize;
      adaptive_quickselect(&mut r, k);
      let answer = r[k];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[k], answer);
    }
  }
}
