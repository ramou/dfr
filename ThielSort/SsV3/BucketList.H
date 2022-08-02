/* $Header:  $ */ 

#ifndef BucketListH
#define BucketListH

/********************************************************************
 *
 *                            BucketList.H
 *
 * Maintain an unordered list of buckets, stored as an array of
 * shorts. The template parameter N_Buckets is the maximum number of
 * buckets in the list. If template parameter Test_Inclusion is true,
 * tag each item added to the list so we can test if a value bkt has
 * been included with GetTag().
 ********************************************************************/

// BucketList                                                    BucketList
template<int N_Buckets=256, bool Test_Inclusion = true>
class BucketList {
  int NBkts;  // Number of buckets in this list.
  int List[N_Buckets];
  // TagOfBkt[b] is 1 iff b is in List, 0 otherwise.
  UByte TagOfBkt[N_Buckets * Test_Inclusion];  
public:
  BucketList() { Reset(); }
  // Replace the bucket in List[ndx] with the last one and shorten
  // the list.
  void Delete(int ndx) { List[ndx] = List[--NBkts]; }
  // Add bkt to the list.
  void Append(short bkt) {
    List[NBkts++] = bkt;
    if (Test_Inclusion) TagOfBkt[bkt] = 1;  // Needed????
  }
  // Add bkt to the list if inc == -1, don't add it if inc ==
  // 0. Behaviour for other values o inc is undefined!!
  void Append(short bkt, int inc) {
    List[NBkts] = bkt;
    NBkts -= inc;
    if (Test_Inclusion) TagOfBkt[bkt] = 1;  // Needed????
  }
  void Reset() {
    NBkts = 0;
    if (Test_Inclusion)
      for (int i = N_Buckets - 1; i >= 0; i--) TagOfBkt[i] = 0;
  }
  int GetN() const { return NBkts; }
  int *GetList() { return List; }
  int GetItem(int ndx) const { return List[ndx]; }
  int GetTag(int ndx) const {
    if (!Test_Inclusion) MyAbort("Invalid call to BucketList::GetTag()");
    return TagOfBkt[ndx];
  }
};  // BucketList

#endif
