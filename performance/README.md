I'd like to give an overview of performance improvements. They're notable. I show the performance in terms of log<sub>10</sub> because otherwise you don't see how it scales so well. Raw data is available. Time is in microseconds. `DFR` is just using `perform` from `make perform`, `RADULS` is using `./raduls -rec_size8 -key_size8 -n_threads1` with the appropriate `-n_recs#` and `Ska Sort` is getting called with:
```C++
radix_sort(values, values + targetLength, buffer,detail::IdentityFunctor());
```
as I couldn't figure out how to use it's initial clone to trivially give what I wanted to know. 

`RADULS` supports multithread, but we're only comparing single-thread here, on the grounds that if you're sorting things all the time, you don't get much ground using a thread that could be sorting something else without overhead... but it's worth noting in case you wanted that absolute speed bonus, because it is there!.

## Titan
I tested each value once and just grabbed whatever. I hope to do proper testing eventually, but these numbers are in-range.
![Titan](https://github.com/ramou/dfr/blob/2acc01ccc589e5e467c4f6d579968595d1819a71/performance/titan.png)

[Source Data](https://github.com/ramou/dfr/blob/2acc01ccc589e5e467c4f6d579968595d1819a71/performance/titan.txt)
