/* $Header:  $ */ 

#ifndef MakeRandomH
#define MakeRandomH

#include <bits/stdc++.h>

/********************************************************************
 *
 *                            MakeRandom.h
 *
 * Generate n random numbers of type UINT (unsigned int or unsigned
 * long) and store them into data. 
 *
 * If true_random is true, create a random seed, otherwise use
 * std::mt19937::default_seed. 
 *
 * If n_bits_to_sort <= 0, all bits are used. If n_bits_to_sort > 0,
 * we use a mask to zero bits above the low order n_bits_to_sort.
 *
 * If n_bits_to_sort > 0 and print_mask_msg == true, print the
 * n_bits_to_sort value to cout.
 *
 * Return the unsigned seed that was used.
 ********************************************************************/

template<typename UINT>
unsigned MakeRandomVec(UINT *data, long N, bool true_random = true,
		   int n_bits_to_sort = 0, bool print_mask_msg = true) {
  UINT max = std::numeric_limits<UINT>::max();
  std::random_device rd;
  auto seed = rd();
  if (!true_random) seed = std::mt19937::default_seed;
  std::mt19937 generator(seed);
  std::uniform_int_distribution<UINT> distribution(0,max);
  if (n_bits_to_sort <= 0) {
    for(unsigned i = 0; i<N;i++)data[i]=distribution(generator);
  } else {
    if (print_mask_msg) 
      std::cout << "Using N_BITS_TO_SORT = " << n_bits_to_sort << '\n';
    UINT mask = 0;
    mask = (~mask) >> (sizeof(UINT) * 8 - n_bits_to_sort);
    for(unsigned i = 0; i<N;i++)data[i]=distribution(generator) & mask;
  }
  return seed;
}
template<typename UINT>
unsigned MakeNormalRandomVec(UINT *data, long N, double avg, double dev,
			     int seed_val = 0) {
  UINT max = std::numeric_limits<UINT>::max();
  std::random_device rd;
#ifdef DEFAULT_SEED
  auto seed = std::mt19937::default_seed;
#else
  auto seed = (seed_val == 0) ? rd() : seed_val;
#endif
  std::mt19937 generator(seed);
  std::normal_distribution<> distribution{avg, dev};
  for(unsigned i = 0; i<N;++i) data[i]=distribution(generator);
  return seed;
}
template<typename UINT>
void MakeSeededRandomVec(UINT *data, long N, unsigned use_seed,
                   int n_bits_to_sort = 0, bool print_mask_msg = true) {
  UINT max = std::numeric_limits<UINT>::max();
  std::mt19937 generator(use_seed);
  std::uniform_int_distribution<UINT> distribution(0,max);
  if (n_bits_to_sort <= 0) {
    for(unsigned i = 0; i<N;i++)data[i]=distribution(generator);
  } else {
    if (print_mask_msg)
      std::cout << "Using N_BITS_TO_SORT = " << n_bits_to_sort << '\n';
    UINT mask = 0;
    mask = (~mask) >> (sizeof(UINT) * 8 - n_bits_to_sort);
    for(unsigned i = 0; i<N;i++)data[i]=distribution(generator) & mask;
  }
}


  
#endif
