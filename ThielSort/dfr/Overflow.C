/******************************************************************
 *
 *                       Overflow.C
 *
 * This file is included within class Overflow and contains private
 * functions. The names with a "_imp" suffix are called by public
 * transfer functions.
 *
 * When a bucket overflows during the deal this class allocates a
 * block of the required size for the data to be stored into. At the
 * end of a pass, it initializes the bucket overflow areas and moves
 * the overflow data into them.
 *
 * Unlike the source and final destination bucket buffers, which
 * belong to the real world, the overflow buffers are used only during
 * a pass.  Their addrs are only used to initialize target vectors.
 * There are 3 overflow buffers:
 *
 * - preface overflow buffer : holds undelimited overflows streams.
 *   Filled during the pass over the source and emptied when the
 *   bucket overflow areas are defined at the end of the pass. A
 *   single buffer local to Overflow sufficies for all passes because
 *   it is emptied before the next pass starts. It may need to be
 *   enlarged at the start of a pass. It has a length set at the start
 *   of the pass, but may run out of space, so GetPrefaceBlock() may
 *   return nullptr.
 *
 * - primary overflow buffer : We use the source space for this. Like
 *   the preface, it is emptied before the end of the pass but the
 *   space is managed by the deal caller. Since we use the source
 *   vector (as it is emptied) it can never run out of space, so
 *   GetPrimaryBlock() always returns a pointer to a block of the
 *   requested size.
 *
 * - merged overflow buffer : the vectors delimiting overflow areas
 *   for individual buckets specify parts of this buffer. Thus it
 *   holds part of the input stream for the next pass. Like the
 *   preface area, it must be extensible (between passes). And like
 *   the preface area, its data will be comsumed before the next pass
 *   needs to start filling it, so a single copy local to Overflow
 *   suffices.
 *******************************************************************/

// OverflowInit_imp                                       OverflowInit_imp

void OverflowInit_imp(Elem_T *over_flow_safe, int safety_margin) {
  SafeAreaOrigin = over_flow_safe;

  // Allocate a preface space even on first pass to protect against
  // pathological case. When RtPerMu > 1, if everything overflows due
  // to estimates of 0, we may run out of space in the primary area if
  // we go there first and allocate blocks of RtPerMu for leftovers.
  // Unlikely in real cases, but happens when n is small (e.g.,
  // 200). Possible in principle for any n I think.
  PrefOfMem.GetAdjustedAdr(&PrefaceOrigin, MergedUsed + safety_margin);
  CurAreaNext = PrefaceOrigin;
  CurAreaEnd = PrefaceOrigin + MergedUsed;
  SafeAreaEnd = SafeAreaOrigin;
  UsingPreface = true;
  
}  // OverflowInit_imp
// Move2Bkts                                                       Move2Bkts

template<typename DBA>
void Move2Bkts(const Elem_T *area_vec, const Elem_T *area_end, int align_size) {
#if 0
  std::cout << "Move2Bkts moving " << area_end - area_vec << std::endl;
#endif
  while (area_vec < area_end) {
    const Elem_T *start = area_vec;
    int t = DBA::GetDealB(area_vec);
    for (area_vec += align_size;
	 (area_vec < area_end) && (t == DBA::GetDealB(area_vec));
	 area_vec += align_size) ;
    int n = area_vec - start;
#ifdef OCD_TESTING
    if (OFBktV[t] == nullptr) {  // Deal error!
      std::cout << "t=" << t << ", n=" << n << ", start = " << start
		<< ", OFBktV[t] = " << OFBktV[t] << "\n";
      ErrTag("Bad value in overflow", true);
    }
#endif
    MUCopy(OFBktV[t], start, n);
    OFBktV[t] += n;
  }
}  // Move2Bkts

// Max amount of padding we might need in target vector if using
// compression technique. 
static constexpr int MaxPadding = (2 * NBuckets * RtPerMu);

// AllocOfBkts                                              AllocOfBkts

void AllocOfBkts(int *of_bkt_list, int n_ofbkts,
		 Elem_T** bkt_addr, const CountType *bkt_n,
		 const int *of_counts, int n_ofe, int align_size) {
#if 0
  std::cout << "AllocOfBkts_imp (RecursiveLevel = " << RecursiveLevel <<
    "): N overflowing bkts = " << n_ofbkts <<
    ", n_ofe = " << n_ofe;
  std::cout << " bkt_addr[0] = " << bkt_addr[0] << std::endl;
#endif
  if (UsingPreface) {  // Preface is still current OF area
    PrefaceEnd = CurAreaNext;
    SafeAreaEnd = SafeAreaOrigin;
  } else {
    SafeAreaEnd = CurAreaNext;
  }
  
  for (int i = 0; i < N_Buckets; i++) OFBktV[i] = nullptr;

  // If overflow > compressThreshold, we relocate the pbuckets and place
  // most chunks in the target buffer.
  CountType compressThreshold = SafeN * OF_COMPRESS_THRESHOLD;
  if (n_ofe > compressThreshold) {  // Fit OF into primary space.
    OfChunksInTarget(of_bkt_list, n_ofbkts, bkt_addr, bkt_n,
		     of_counts, n_ofe, bkt_addr[0]);
  } else {
    OfChunksInMergedOf(of_bkt_list, n_ofbkts, of_counts, n_ofe, align_size);
  }
}

/**
 * Relocate the bkts so almost all the data will fit in the primary
 * area instead of the OF area. Set the pbkt addrs and move the pbkts
 * to their final location. Set the ofbkt addrs, but don't do the
 * redeal here.
 **/
void OfChunksInTarget(int *of_bkt_list, int n_ofbkts,
		      Elem_T** bkt_addr, const CountType *bkt_n,
		      const int *of_counts, int n_ofe,
		      Elem_T* target_bfr) {
  ReorderBkts.InitBktVec(bkt_n, of_counts);
  CountType nOF = ReorderBkts.CreateAreaMaps();  
  MergedOfMem.GetAdjustedAdr(&MergedOrigin, nOF);  // Wrong!! Too big!!  
  ReorderBkts.SetOfChunkAddrs(bkt_addr, OFBktV, target_bfr, MergedOrigin);
}

/**
 * Copy n recs from src_p to dest_p. Use an efficient move if possible.
 * n must be x * RtPerMu.
 **/
void MUCopy(Elem_T *dest_p, const Elem_T *src_p, CountType n) {
  if ((RtPerMu > 1) && MoverP) {  // Invoke flush mover
    Mover_T::CopyMus(const_cast<Elem_T *>(src_p), dest_p, n / RtPerMu);
  } else {
    std::memcpy(dest_p, src_p, n*sizeof(Elem_T));
  }
}
void OfChunksInMergedOf(const int *of_bkt_list, int n_ofbkts,
			const int *of_counts, int n_ofe, int align_size) {
  const int aM1 = align_size - 1;
  // Allocate merged of area. Allow for max alignment adjustment possible.
  //std::cout << "OfChunksInMergedOf requests " << n_ofe << " + slop\n";
  if (n_ofe)  // Last pass if n_ofe is 0.
    MergedOfMem.GetAdjustedAdr(&MergedOrigin, n_ofe + n_ofbkts * aM1);
  Elem_T *OfP = MergedOrigin;
  // Define the target space for each bucket
  for (int i = 0; i < n_ofbkts; i++) {
    int b = of_bkt_list[i];
    OFBktV[b] = OfP;
    OfP += (of_counts[b] + aM1) & ~ aM1;
  }
  MergedUsed = OfP - MergedOrigin;

}

/**
 * Define (but don't fill) the primary and overflow chunks in *pool_p.
 **/
template<typename Pool_T>
static void CreateDpChunks(Elem_T **bkt_addr, const int* bkt_n,
			   Elem_T **ofbkt_addr, const int* ofbkt_n,
			   Pool_T *pool_p) {
  for (int b = 0; b < NBuckets; b++) {
    // Don't need if unless other code checks addr==nullptr for empty.
    // ONce this is working, simplify it.
    if (ofbkt_n[b] > 0) {  // Full bucket + overflow
      pool_p->DefineChunk(b, bkt_addr[b], bkt_n[b],
			  ofbkt_addr[b], ofbkt_n[b]);
    } else {  // OK if len == 0.
      pool_p->DefineChunk(b, bkt_addr[b], bkt_n[b],
			  nullptr, 0);
    }
  }
}
