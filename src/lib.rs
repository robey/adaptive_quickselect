extern crate rand;

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

/// partition by the element at `k`, and return the new position of that element.
pub fn hoare_partition<A: PartialOrd>(r: &mut [A], k: usize) -> usize {
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




/// select the `k`th item from a hypothetically-sorted version of `r`.
/// for example, the median of `r` is `adaptive_quickselect(r, r.len() / 2)`.
/// `r` will be modified by partially sorting (partitioning) around various
/// pivots until the `k`th item is in the right place. no guarantee is made
/// about any other position.
pub fn adaptive_quickselect<A: PartialOrd + std::fmt::Debug>(r: &mut [A], k: usize) {
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
    hoare_partition(r, k)
  } else {
    0
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

  fn make_random_sequence(generator: &mut StdRng, n: usize) -> Vec<f64> {
    let mut rv = Vec::<f64>::with_capacity(n);
    for _ in 0..n {
      rv.push(generator.next_f64());
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
}
