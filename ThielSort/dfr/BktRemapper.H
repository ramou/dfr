/* $Header:  $ */ 

#ifndef BktRemapperH
#define BktRemapperH

/********************************************************************
 *
 *                            BktRemapper.H
 *
 *
 ********************************************************************/

// BktRemapper                                                  BktRemapper
template <typename Elem_T>
struct BktRemapper {
  // Original addr of pbkt. nullptr if empty. Is set to nullptr when
  // the bkt is successfully moved to its remapped location.
  Elem_T* OldPbAdr;
  CountType PBktOffset;  // nullptr or final addr of the pbkt.
  CountType OfBktOffset;  // nullptr or final addr of the ofbkt.
  CountType PbN;  // pbucket size (including padding).
  CountType OfN;  // ofbucket size (including padding).
  int Status;  // 0:neither in OF, 1: pbkt in OF, 2: ofbkt in OF
  bool MStatus;  // false until pbkt is in remapped location.

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  //BktRemapper();
  //~BktRemapper();

  // Init object to empty.
  void InitBktClass() {
    PBktOffset = OfBktOffset = PbN = OfN = 0; Status = 0; OldPbAdr = 0;
    MStatus = false;
  }
  void InitOldAddr(Elem_T* old_p) { OldPbAdr = old_p; }
  /**
   * Return the addr of the ofbucket (OfBkt) or pbucket (PBkt) by
   * adding the bucket offset to the appropriate buffer addr. Return
   * nullptr if the bucket is empty.
   **/
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  /**
   * Return -1 if the space b_p:b_end is entirely before this pbkt, 0
   * if it overlaps this pbkt, and 1 if it is entirely after this
   * pbkt.
   **/
  int OverlapStatus(Elem_T* p, Elem_T* p_end) {
    if (p_end < OldPbAdr) return -1;
    if (p >= (OldPbAdr + PbN)) return 1;
    return 0;
  }
  Elem_T* OfBkt(Elem_T* p_bfr_p, Elem_T* of_bfr_p) {
    if (OfN) return ((Status == 2) ? of_bfr_p : p_bfr_p) + OfBktOffset;
    return nullptr;
  }
  Elem_T* PBkt(Elem_T* p_bfr_p, Elem_T* of_bfr_p) {
    if (PbN) return ((Status == 1) ? of_bfr_p : p_bfr_p) + PBktOffset;
    MStatus = true;  
    return nullptr;
  }
  // Record p_offset as offset of the pbkt and of_offset as the
  // offset for the ofbkt. status_val records the areas holding the
  // chunks. See Status doc.  The reference parameters are both
  // inc'ed. Note that they may reference the same value!!
  void SetOffsets(CountType &p_offset, CountType &of_offset,
		  int status_val) {
    PBktOffset = p_offset; p_offset += PbN;
    OfBktOffset = of_offset; of_offset += OfN;
    Status += status_val;
  }
  /**
   * Handle move to the OF area. 
   **/
  void EasyMove(Elem_T* dest_p) {
    if (!MStatus && (Status ==  1)) {
      DoMove(dest_p);  // Move to OF area is safe
    }
  }
  void LeftMove(Elem_T* dest_p) {
    if (!MStatus && (Status !=  1)) {
      if (dest_p < OldPbAdr) DoMove(dest_p);  // Move left.
    }
  }
  void RightMove(Elem_T* dest_p) {
    if (!MStatus && (Status !=  1)) {
      if (dest_p > OldPbAdr) DoMove(dest_p);  // Move right.
    }
  }
  void DoMove(Elem_T* dest_p) {
    // Skip fast move for now. 
    std::memmove(dest_p, OldPbAdr, PbN*sizeof(Elem_T));
    MStatus = true;
  }
};  // BktRemapper

#endif
