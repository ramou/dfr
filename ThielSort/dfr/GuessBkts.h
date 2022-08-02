/**
 *                         GuessBkts.h
 *
 **/

/**
 * Divide buffer_len ETypes as evenly as possible among n_bkts buckets
 * (while making each bucket start at x*bkt_alignment from 0). Store
 * the resulting bucket lens in bkt_len_v, and sum(blk lens) in
 * bkt_len_v[n_bkts]. 
 *
 * NOTE: If buffer_len is smallish compared to n_bkts, bkt_alignment
 * should be made small enough to allow as many buckets as possible to
 * have nonzero lengths. I can't enforce that here because the calling
 * program may assume alignments, but having a few buckets with
 * bkt_alignment space and most with 0 space probably isn't a good
 * choice. It will give correct results, but many buckets will overflow
 * immediately.
 **/
NO_INLINE
void GuessBkts(int *bkt_len_v, int buffer_len, int bkt_alignment, int n_bkts) {
  int bktLen = buffer_len / n_bkts;  // tentative
  bktLen -= bktLen % bkt_alignment;  // round down to alignment
  int waste = buffer_len - bktLen * n_bkts;
  // Use up as much waste as possible by extending first few buckets
  // by bkt_alignment.
  int nCLAvail = waste / bkt_alignment;
  for (int b = 0; b < nCLAvail; b++) {
    bkt_len_v[b] = bktLen + bkt_alignment;
  }
  for (int b = nCLAvail; b < n_bkts; b++) {
    bkt_len_v[b] = bktLen;
  }
  bkt_len_v[n_bkts] = nCLAvail * (bktLen + bkt_alignment) +
    (n_bkts - nCLAvail) * bktLen;
}

