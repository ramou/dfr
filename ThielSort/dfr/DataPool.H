/* $Header:  $ */ 

#ifndef DataPoolH
#define DataPoolH

#include "AlignedBlock.H"

/********************************************************************
 *
 *                            DataPool.H
 *
 * We refer to the bucket structure that holds the data during
 * intermediate passes of the sort (i.e., the vectors delimiting
 * individual buckets and individual overflow areas), as a
 * DataPool. It is easy to extend this idea to include the initial
 * source vector (by making a data pool with a single "bucket" and no
 * overflow area) and the final result (when the buckets may delimit
 * parts of the final vector). This provides a single simple interface
 * for handling the data both as a source and a destination for all
 * passes. The interleaving of primary bucket areas and overflow
 * area is messy, so it is hidden here.
 *
 * During the sort we swap the data back and forth between a pair of
 * DataPools.
 *
 * Note that the bucket vectors point into 1 data area, but that data
 * area isn't really part of the data pool as its address is only used
 * to initialize bucket vectors at the start of a pass.
 *
 ********************************************************************/

// DataPool                                                       DataPool
template<typename SRec_T>
class DataPoolClass<SRec_T, CHUNK_DATAPOOL> {

  typedef SRec_T ElemT;
  typedef SRec_T BucketT;
  
 public:

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  DataPoolClass() {  }
  ~DataPoolClass() { ReleaseBB(); }

  /**
   * (Re)allocate or initialize the bucket buffer to hold buffer_len
   * ElemTs, and return the addr of the buffer. If bkt_buffer =
   * nullptr, allocate a buffer of buffer_len ElemTs (or reuse the
   * previously allocated one). Otherwise, just use the one
   * given. buffer_len satisfies any alignment requirements.
   **/  
  ElemT *InitBuffer(ElemT *bkt_buffer, long buffer_len,
		    int n_buckets) {
    NChunks = NBuckets = n_buckets;
    if (bkt_buffer != nullptr) {
      BktBuffer = bkt_buffer;
    } else {
      if (PreferNewBuffer) BucketBuffer.PreferCleanSpace();
      BucketBuffer.GetAdjustedAdr(&BktBuffer, buffer_len);
      BBOwner = true;
    }
    BBLen = buffer_len;
    return BktBuffer;
  }
  /**
   * Version for use when segmenting the sort. Allocate a buffer of
   * size buffer_len as the initial buffer. It will be increased if
   * needed. The above version still needs to be called for each
   * segment.  This allows us to allocate a buffer long enough to hold
   * the biggest segment instead of having to keep reallocating longer
   * buffers as we go through the segments.
   **/
  void InitBuffer(long buffer_len) {
    BBLen = buffer_len;
    if (PreferNewBuffer) BucketBuffer.PreferCleanSpace();
    BucketBuffer.GetAdjustedAdr(&BktBuffer, buffer_len);
    BucketBuffer.MarkFree();
  }

  
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.
  
  // Return the current buffer.
  ElemT *GetBktBuffer() { return BktBuffer; }
  //// Deal the DataPool

  // Do successive passes. Source was created in this dp by the
  // previous pass, or wrapped here with WrapSourceVec().
  void DealFromDP(DealerBase<BucketT> *dealer_p) {
    //CountType RecsToDealer = 0;   // DEBUG
    for (int b = 0; b < NChunks; b++) {
      if (PrimaryChunkV[b].Len > 0) {
	//RecsToDealer += PrimaryChunkV[b].Len;  // DEBUG
	dealer_p->ProcessBlk(PrimaryChunkV[b].Addr, PrimaryChunkV[b].Len);
      }
      if (OfChunkV[b].Len > 0) {
	//RecsToDealer += OfChunkV[b].Len;  // DEBUG
	dealer_p->ProcessBlk(OfChunkV[b].Addr, OfChunkV[b].Len);
      }
    }
    //std::cout << RecsToDealer << " recs passed to dealer\n";  // DEBUG
    ReleaseBB();
  }
  void ReleaseBB() {
    if (BBOwner) { 
      BucketBuffer.MarkFree();
      BBOwner = false;
    }
  }    
  /**
   * Wrap the source vector src_v : src_v + src_n in this DataPool for
   * first pass.
   **/
  void WrapSourceVec(ElemT *src_v, long src_n) {
    PrimaryChunkV[0].Store(src_v, src_n);
    OfChunkV[0].Store(nullptr, 0);
    NChunks = 1;
    BktBuffer = src_v;  // Safe overflow area
  }
  /**
   * Copy the data to a contiguous target vector. dest_vec is assumed
   * to be large enough to hold the data.
   *
   * NOTE: Can we do better with simd? But alignment is a complication...
   **/
  void CopyDpToVec(ElemT *dest_vec) {
    for (int b = 0; b < NChunks; b++) {
      if (PrimaryChunkV[b].Len > 0) {
 	std::memcpy(dest_vec, PrimaryChunkV[b].Addr,
		    PrimaryChunkV[b].Len*sizeof(ElemT));
	dest_vec += PrimaryChunkV[b].Len;
      }
      if (OfChunkV[b].Len > 0) {
	std::memcpy(dest_vec, OfChunkV[b].Addr,
		    OfChunkV[b].Len*sizeof(ElemT));
	dest_vec += OfChunkV[b].Len;
      }
    }
    ReleaseBB();
  }
  /**
   * For debugging: compute and print the sum of the primary and OF
   * chunks (separately).
   **/
  void CountRecs(CountType expected_n, const char *label) {
    CountType n = 0;
    CountType nOF = 0;
    for (int b = 0; b < NChunks; b++) {
      n += PrimaryChunkV[b].Len;
      nOF += OfChunkV[b].Len;
    }
    std::cout << label;
    std::cout << n << " recs in primary and " << nOF <<
      " recs in OF chunks (total = " << n + nOF << ").";
    if ((n + nOF) == expected_n) std::cout << " OK.\n";
    else std::cout << "***Error! Should be " << expected_n << '\n';
  }
  /**
   * Use CopySort() to combine copy after odd number of deal passes to
   * dest_v with final insertion sort.
   **/
  void InsertionSortWithCopy
  (InsertionSortPass<ElemT> *is_obj, ElemT top_mask, ElemT *dest_v) {
    int b = 0;
    while (is_obj->CopySort(PrimaryChunkV[b].Addr, PrimaryChunkV[b].Len,
			    OfChunkV[b].Addr, OfChunkV[b].Len,
			    top_mask, dest_v)) b++;
    // Have copied / sorted first non empty bucket.
    for (b++; b < NBuckets; b++) {
      is_obj->ExtendCopySort(PrimaryChunkV[b].Addr, PrimaryChunkV[b].Len,
			     OfChunkV[b].Addr, OfChunkV[b].Len);
    }
    ReleaseBB();
  }
  
// ConfirmTargetBkts                                      ConfirmTargetBkts
  /**
   * Check datapool and ensure that chunks sum to expected_n and that
   * each chunk has the same Deal Byte (from DBA). Always writes a
   * message.
   **/
  template<typename DBA> void ConfirmTargetBkts(long expected_n) {
    const char* ChnkType[] = { "For primary bucket", "For Overflow bucket" };
    int nErrs = 0;
    int dealB = 0;
    int nPrimChunks = 0;
    int ChnkCnts[] = { 0, 0};
    using namespace std;
    for (int b = 0; b < NChunks; b++) {
      ElemT *p = PrimaryChunkV[b].Addr;
      long n = PrimaryChunkV[b].Len;
      for (int chnk = 0; chnk < 2; chnk++) {  // chnk = primary, overflow
	if (n > 0) {
	  ChnkCnts[chnk] += n;
	  int LclNErrs = 0;
	  for (int j = 0; j < n; j++) {
	    int curB = DBA::GetDealB(p + j);
	    if (b != curB) {
	      LclNErrs++;
	      if (LclNErrs <= 8) BadDealByte(b, ChnkType[chnk], j, curB);
	    }
	  }
	  if (LclNErrs > 0) {
	    cout << ChnkType[chnk];
	    InBrackets(b);
	    cout <<  " had " << LclNErrs << " bad elements\n";
	    nErrs += LclNErrs;
	  }
	}
	n = OfChunkV[b].Len;
	p = OfChunkV[b].Addr;
      }
    }
    int TotalN = ChnkCnts[0] + ChnkCnts[1];
    cout << ChnkCnts[0] << " items in buckets and " << ChnkCnts[1] <<
      " items in overflow areas.\n";
    if ((TotalN != expected_n) || (nErrs != 0)) {
      if (TotalN != expected_n)
	cout << "Wrong count: expected " << expected_n << " found " <<
	  TotalN << "\n";
      else cout << "Counts are consistent. ";
      if (nErrs == 0) cout << "Chunks are consistent.\n";
      else {
	cout << " consistency errors within chunks found.\n";
      }
    } 
  }
  /**
   * Store values for primary and overflow chunks for bucket bkt_id.
   **/
  void DefineChunk(int bkt_id, BucketT *primry_p, int primry_n,
		   BucketT *of_p, int of_n) {
    PrimaryChunkV[bkt_id].Store(primry_p, primry_n);
    OfChunkV[bkt_id].Store(of_p, of_n);
  }
private:
  /**
   * The data in a particular bucket sits in either 1 or 2 chunks,
   * depending on whether an overflow area is used. We use a vector of
   * chunks (PrimaryChunkV) to describe the location of the data in
   * the main bucket buffer, and a second, parallel vector (OfChunkV)
   * to describe any data that has been placed in the overflow area. We
   * have 3 cases:
   *
   * - bucket b is empty. PrimaryChunkV[b].Len is
   *   0. PrimaryChunkV[b].Addr and OfChunkV[b] are undefined.
   *
   * - bucket b has data but doesn't need the overflow area.
   *   PrimaryChunkV[b].Addr points to the start of the main data
   *   block, PrimaryChunkV[b].Len is the number of items in the
   *   bucket.  OfChunkV[b].Len is 0, OfChunkV[b].Addr is undefined.
   *
   * - bucket b has data and is overflowing.  PrimaryChunkV[b].Addr
   *   points to the start of the main data block,
   *   PrimaryChunkV[b].Len is the number of items in the main data
   *   block, OfChunkV[b].Addr points to the start of the overflow
   *   area for bucket b, OfChunkV[b].len is the number of items in
   *   the overflow area for bucket b.
   **/
  struct ChunkT {
    BucketT *Addr;
    long Len;
    void Store(BucketT *p, long n) {
#ifdef OCD_TESTING
      if ((p == nullptr) & (n != 0)) ErrTag("Bad store in ChunkT");
#endif      
      Addr = p; Len = n;
    }
  };
  ChunkT PrimaryChunkV[MaxNBuckets];
  ChunkT OfChunkV[MaxNBuckets];
  
  BucketT *BktBuffer = nullptr;    // The start of the main bucket buffer
  long BBLen = 0;         // Actual length of BktBuffer
  bool BBOwner = false;
  int NBuckets;
  // N bucket-pair chunks in current datapool. Normally == NBuckets,
  // but will be 1 if this DataPool is wraping a source vector.
  int NChunks;  
  
};  // DataPool

#endif
