extern crate rand;

use std::fmt::Debug;
// #[macro_use]
// extern crate lazy_static;

#[inline]
pub fn median3<A: PartialOrd>(r: &mut [A], a: usize, b: usize, c: usize) {
  if r[b] < r[a] {
    if r[b] < r[c] {
      if r[c] < r[a] {
        // b < c < a
        r.swap(b, c);
      } else {
        // b < a <= c
        r.swap(b, a);
      }
    }
    // else: c <= b < a
  } else if r[c] < r[b] {
    if r[c] < r[a] {
      // c < a <= b
      r.swap(b, a);
    } else {
      // a <= c < b
      r.swap(b, c);
    }
  }
  // else: a <= b <= c
  debug_assert!((r[a] <= r[b] && r[b] <= r[c]) || (r[a] >= r[b] && r[b] >= r[c]));
}




/// Returns the index of the median of r[a], r[b], r[c].
#[inline]
pub fn median_index<A: PartialOrd>(r: &[A], a: usize, b: usize, c: usize) -> usize {
  if r[a] > r[c] {
    if r[b] > r[a] { a } else { if r[b] < r[c] { c } else { b } }
  } else {
    if r[b] > r[c] { c } else { if r[b] < r[a] { a } else { b } }
  }
}

/// Tukey's ninther: compute the median of each triplet at a time, then the
/// median of those medians, returning that index.
pub fn ninther<A: PartialOrd>(
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
  // println!("--> expand_right pivot={:?} hi={:?} right={:?}", pivot, hi, right);
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

/// Same as expand_partition_right, but in reverse.
/// [ < < x x x x < < < p ]
///       left    lo    pivot
fn expand_partition_left<A: PartialOrd>(r: &mut [A], pivot: usize, lo: usize, left: usize) -> usize {
  // println!("--> expand_left pivot={:?} lo={:?} left={:?}", pivot, lo, left);
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
  // println!("--> expand {:?}", r);
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
pub fn partition_hoare<A: PartialOrd>(r: &mut [A], k: usize) -> usize {
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

/// partition by finding a pivot that approximates the median.
pub fn partition_ninthers<A: PartialOrd + Debug>(r: &mut [A]) -> usize {
  // println!("--> Ninthers: {:?}", r);
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
  // println!("    now: {:?}", r);
  // println!("    pivot={:?} lo={:?} hi={:?}", pivot, lo, hi);

  adaptive_quickselect(&mut r[lo..hi], pivot);
  let x = expand_partition(r, lo, lo + pivot, hi);
  // println!("    now {:?}", r);
  // println!("    pivot={:?}", x);
  x
}







/// select the `k`th item from a hypothetically-sorted version of `r`.
/// for example, the median of `r` is `adaptive_quickselect(r, r.len() / 2)`.
/// `r` will be modified by partially sorting (partitioning) around various
/// pivots until the `k`th item is in the right place. no guarantee is made
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



#[cfg(test)]
mod tests {
  use adaptive_quickselect;
  use median3;
  use rand::{Rng, SeedableRng, StdRng};

  fn make_random_sequence(generator: &mut StdRng, n: usize) -> Vec<u16> {
    let mut rv = Vec::<u16>::with_capacity(n);
    for _ in 0..n {
      rv.push((generator.next_u32() & 0xffff) as u16);
    }
    rv
  }

  #[test]
  fn test_median3() {
    let mut t1 = [ 3, 4, 5 ];
    median3(&mut t1, 0, 1, 2);
    assert_eq!(t1[1], 4);

    let mut t2 = [ 5, 4, 3 ];
    median3(&mut t2, 0, 1, 2);
    assert_eq!(t2[1], 4);

    let mut t3 = [ 5, 3, 4 ];
    median3(&mut t3, 0, 1, 2);
    assert_eq!(t3[1], 4);

    let mut t4 = [ 4, 5, 3 ];
    median3(&mut t4, 0, 1, 2);
    assert_eq!(t4[1], 4);
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
  fn test_ninthers() {
    let mut generator = StdRng::from_seed(&[ 1337 ]);

    for _ in 0..5000 {
      let mut r = make_random_sequence(&mut generator, 1024);
      adaptive_quickselect(&mut r, 512);
      let answer = r[512];
      r.sort_by(|a, b| a.partial_cmp(b).unwrap());
      assert_eq!(r[512], answer);
    }
  }
}
