/* $Header:  $ */ 

#ifndef DataInfoH
#define DataInfoH

/********************************************************************
 *
 *                            DataInfo.H
 *
 * A class to collect info about the data being dealt. Currently
 * LiveBits and counts (or nothing) can be collected. The digit to be
 * counted depends on the Count_Byte template parameter in the DBA
 * template argument passed to Collect<>(). Collect<>() does both,
 * controlled by a compile time test of the Action_Options template
 * parameter. When no actions are turned on, the compiler should
 * delete the code completely, as well as the parameter evaluation for
 * parameters no longer used.
 * 
 * See the templated "using" declarations at the bottom for standard
 * template names for the various collect cases.
 ********************************************************************/

enum DataInfoOptions
  { DiNoOptions = 0, DiLbOption = 1, DiCountOption = 2 };

template<typename Elem_T, int Action_Options>class DataInfo;

// DataInfo                                                       DataInfo
template<typename Elem_T, int Action_Options>
class DataInfo {
  void ResetLiveBit(Elem_T init_bitmask) {
    if (Action_Options & DiLbOption) {
      LiveBits = 0;
      BitMask = init_bitmask;
    }
  }
  void ResetCountVec(CountType *init_counts, int n_buckets) {
    if (Action_Options & DiCountOption) {
      Counts = init_counts;
      for (int i = 0; i < n_buckets; i++) Counts[i] = 0;
    }
  }
protected:
  CountType *Counts = nullptr;  // Accumulate count here if counting.
  Elem_T LiveBits = 0;
  Elem_T BitMask = 0;
  
 public:
  static const int DataInfoOptions = Action_Options;
  
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  DataInfo(Elem_T init_bitmask = 0, CountType *init_counts = nullptr) :
    Counts(init_counts), BitMask(init_bitmask) {}
  //~DataInfo();

  /**
   * Since we may need to reuse a DataInfo object for several subsorts
   * when splitting a long sort into more manageable pieces (e.g.,
   * SdDivertAndConquor.H), we may need to reset the BitMask and / or
   * count vector, as required by the Action_Options.
   **/
  void ResetDataInfo(Elem_T init_bitmask, int n_buckets,
		     CountType *init_counts = nullptr) {
    ResetLiveBit(init_bitmask);
    ResetCountVec(init_counts, n_buckets);
  }
  
  /**
   * Since DataInfo<> is a base class for the Dealer, and the Dealer
   * is a strictly local object, we use an auxillary DataInfo<> object
   * (of the same type) to pass information into and out of the base
   * class, so neither the Dealer nor the DataPoolPair<> object need
   * to see the details. The auxillary object is created in the Sd
   * object and passed through to the base class using the following
   * functions.
   **/
  
  // Init all values from src before starting. Note that the count
  // vector is zeroed here if appropriate.
  void InitDataInfo(const DataInfo &src, int n_buckets) {
    // Used to copy LiveBits too. Now just set to 0.
    ResetLiveBit(src.BitMask);
    ResetCountVec(src.Counts, n_buckets);
  }
  // Extract computed LiveBits from src. Rest not relevant
  // since Counts vector addr is passed in.
  void ExtractInfo(const DataInfo &src) {
    if (Action_Options & DiLbOption) LiveBits = src.LiveBits;
  }
  
  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  Elem_T GetLiveBits() { return LiveBits; }

  template<typename DBA>
  void Collect(const Elem_T &v, Elem_T &live_bits, const Elem_T bit_mask) {
    if (Action_Options & DiLbOption) live_bits |= bit_mask ^ v;
    if (Action_Options & DiCountOption) Counts[DBA::GetCountB(v)]++;
  }
    
};  // DataInfo

// Collect no info during deal.
template <typename Elem_T>
using NoopDataInfoT = DataInfo<Elem_T, DiNoOptions>;
// Collect livebits during deal.
template <typename Elem_T>
using LbDataInfoT = DataInfo<Elem_T, DiLbOption>;
// Collect counts during deal.
template <typename Elem_T>
using CountDataInfoT = DataInfo<Elem_T, DiCountOption>;
// Collect both livebits and counts during deal.
template <typename Elem_T>
using CntLbDataInfoT = DataInfo<Elem_T, DiLbOption + DiCountOption>;

#endif
