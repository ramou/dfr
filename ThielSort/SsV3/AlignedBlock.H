
/* $Header:  $ */ 

#ifndef AlignedBlockH
#define AlignedBlockH

/********************************************************************
 *
 *                            AlignedBlock.H
 *
 * Container to provide a block of memory that has the low order bits
 * of its address set to 0 (i.e., to be aligned on a certain boundry,
 * such as a cache line boundry). The original address is saved and
 * the block automatically deleted when the AlignedBlock object is
 * destroyed. Function GetAdjustedAdr() is used to allocate the
 * memory, adjust it, and return the aligned block.  Each AlignedBlock
 * object can only hold 1 object or array at a time. Additional calls
 * to GetAdjustedAdr() will reuse the same space if it is big enough,
 * or free it and allocate new space.
 *
 * If template parameter Control_Access is true, GetAdjustedAdr() will
 * mark the block as in use (the default is available). MarkFree()
 * marks the block as available, but doesn't delete the allocated
 * block. A call to GetAdjustedAdr() or an attempt to delete the
 * AlignedBlock object while the block is marked busy is a fatal
 * error. This allows the block to be reused safely serially or in a
 * recursion, when it is "known" that the multiple uses won't overlap.
 *
 * Note that Control_Access can't detect access via a saved pointer,
 * even if the block has been marked as free and reallocated. Caveat
 * emptor.
 ********************************************************************/

// AlignedBlock                                               AlignedBlock
template<bool Control_Access = false>
class AlignedBlock {
 private:
  static const bool UseMalloc = true;
  void *Mem = nullptr;
  size_t NBytes = 0;
  bool Inuse = false;
  /**
   * Some machines (e.g., titan) seem to perform better using newly
   * allocated space than reusing old space for a buffer. If
   * AlwaysRealloc is true, we release any current space and
   * reallocate it in AllocSpace(). If requesting the same size, it
   * seems likely that we will get the same block back sometimes
   * (often?), which should (??) be cheaper. But in cases where the
   * required amount varies somewhat, allocating a marginally
   * smaller space may cause a "completely new" block to be
   * created. So if we have AlwaysRealloc true I always (re)allocate
   * a block of size max(needed, NBytes). I don't know what
   * mechanisms are in play here, so this may well be redundant, but
   * I have seen the "prefer new space" scenario affect the cost of
   * a sort on titan often. Something is going on and this seems to
   * help.
   **/
  bool AlwaysRealloc = false;
  const char *WhoAmI = "";  // Empty unless naming contructor called.
  void Die() {
    ErrTag("Attempt to use / delete an inuse block (AlignedBlock)");
    ErrTag(WhoAmI, true);
  }
  void Delete() {
    if (Control_Access & Inuse) Die();
    if (UseMalloc) free(Mem);
    else delete [] (char *)Mem;
  }
  void AllocSpace(size_t needed) {
    if ((needed > NBytes) || AlwaysRealloc) {
      if (Mem) Delete();
      needed = std::max(NBytes, needed);
      if (UseMalloc) Mem = malloc(needed);
      else Mem = new char [needed];
      NBytes = needed;
    }  // else reuse previous space.
  }
  void *Request(long byte_size, int mask = 0x3f) {
    if (Control_Access & Inuse) Die();
    long needed = byte_size + mask + 1;  // ensure enough space after align.
    AllocSpace(needed);
    void *t = (void *)(((long)Mem + mask) & ~mask);
    //cout << "Changing Target from " << Mem << " to " << t << endl;
    if (Control_Access) Inuse = true;
    return t;
  }
 public:

  //////////////////////////////////////////////////////////////////////
  //                   Constructors / Destructor

  // Must create empty! The block allocate and address return are done
  // together in GetAdjustedAdr().
  AlignedBlock() {}
  AlignedBlock(const char *s) : WhoAmI(s) {}  // Who is complaining??
  ~AlignedBlock() { Delete(); }

  //////////////////////////////////////////////////////////////////////
  //                   Other functions
  
  /**
   * Return a block that is byte_size chars long, aligned so that the
   * starting address & mask is 0. Default is aligned on a 64 byte
   * boundry. 
   *
   * NOTE: No constructor or destructor is called for your
   * object(s). You will have to do that yourself. But the memory will
   * be released when the AlignedBlock is destroyed. Do not explicitly
   * release the allocated space.
   **/
  void *GetAdjustedAdr(long byte_size, int mask = 0x3f) {
    return Request(byte_size, mask);
  }
  /**
   * As above, but allocate an aligned vector of type T[n_eles] and
   * store its addr into target_adr.
   **/
  template<class T>
    void GetAdjustedAdr(T **target_adr, long n_eles, int mask = 0x3f) {
    *target_adr = new(Request(n_eles * sizeof(T), mask)) T;
  }
  /**
   * Similar to GetAdjustedAdr() except that the correct pointer type
   * is returned and the default constuctor is called for the object.
   * The object type must be provided as a template parameter. Not for
   * vectors.
   **/  
  template<class T>T* ConstructAlignedObject(int mask = 0x3f) {
    return new(Request(sizeof(T), mask)) T;
  }
  void MarkFree() { Inuse = false; }
  bool IsInuse() { return Inuse; }
#if 0
  // p just to allow T to be identified, since above gives me grief.
  // Have to identify fn with template keyword??
  // *** Deprecated in favor of the above versions.
  template<class T>T* ConstructAlignedObject(T* p, int mask = 0x3f) {
    return new(Request(sizeof(T), mask)) T;
  }
#endif
  // Ensure that next request gets a newly allocated space.
  void PreferCleanSpace() { NBytes = 0; }
};  // AlignedBlock

#endif

