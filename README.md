# adaptive_quickselect

This is a rust implementation of the "adaptive quickselect" algorithm for finding the `k`th-ranked element of an unsorted list in linear time. You could use this to find the median or 25th percentile of a very large array of numbers, for example.


## How?

The partitioning scheme from quicksort can be used to divide-and-conquer a large list of numbers, recursively partitioning whichever half has the desired rank in it. Zeno's race means it should typically converge to 2N time, though it can be N-squared in the worst case. This algorithm is described well here: https://rcoh.me/posts/linear-time-median-finding/

Adaptive quickselect uses some tricks to throw away the worst choices for a partitioning pivot, and guarantee linear time. Andrei Alexandrescu expanded these tricks to cover the case when the desired rank is close to either end of the list as well as the median, and added clever optimizations so that the work of sampling potential pivots can be used to jumpstart the partitioning pass.

  - [Fast Deterministic Selecting](http://erdani.com/research/sea2017.pdf) - a paper describing the algorithms
  - [MedianOfNinthers](https://github.com/andralex/MedianOfNinthers/) - a C++ implementation of `adaptive_quickselect`

All I did was port the code to rust, simplifying it in the process (because rust's syntax is a bit more expressive). I also modified `expand_partition_{left,right}` to use one pass, because careful consideration over several coffees convinced me that splitting the code into two cases was unnecessary.


## Build

```
$ cargo test
```

The tests may take over 10 seconds, depending on your computer. They use a PRNG to create tens of thousands of arrays, select a rank, and then sort the array and verify that the answer was correct.


## Use

```rust
let mut r = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 12, 13, 14, 15, 16, 17, 18, 19, 20, 11 ];
adaptive_quickselect(&mut r, 10);
assert_eq!(r[10], 11);
```


## Author

- Robey Pointer <robeypointer@gmail.com>

This code is licensed under the Apache 2.0 license, found here: http://www.apache.org/licenses/LICENSE-2.0
