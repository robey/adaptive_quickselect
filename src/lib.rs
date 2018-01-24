//use std::mem;

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


#[cfg(test)]
mod tests {
  use median3;

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

}
