/* $Header:  $ */ 

#ifndef ConfigH
#define ConfigH

/**
 * Configuration testing and logging is a mess, given the interactions
 * that can occur. There are two aspects to a configuration: 
 *
 * 1 - integer and boolean constants provided via #defines to control
 *     module behaviour. There are several of them. Most, but not all
 *     of them are relevant for every SD. All of them have defaults,
 *     and some of them have different defaults depending on the
 *     requested SD. After they are all given values, these values are
 *     used to define static const members in namespace ConstConfig.
 *
 * 2 - module selection chooses 1 module from each of the following types:
 *     LADLE_MANAGER
 *     FLUSHER_MANAGER
 *     FLUSH_MOVER
 *     SORT_DESCRIPTOR
 *
 *   The module selection is done by a sort descriptor config class
 *   (SdConfig*).  The CONFIG_PATTERN macro selects a single SdConfig
 *   class with the #if nest below to be added to namespace
 *   ConstConfig and compiled.
 *
 **/

// These includes are needed for almost all config patterns.

// System includes
#include <cassert>
#include <iostream>
#include <bits/stdc++.h>
#include <limits>
#include <algorithm>

// My support utilities not inherently connected with sorting.

#include "OddsAndEnds.h"
#include "AlignedBlock.H"
#include "CmdLineScanner.H"
#include "RandomTool.H"
#include "BucketList.H"
#include "VectorType.H"
#include "PfSupport.H"
#include "ChronoStopWatch.H"

#include "ConstConfig.H"

namespace ConstConfig {

  int RecursiveLevel = 0;
  const CountType DiversionThresholds[DIVERSION_ARRAY_SIZE] =
    { DIVERSION_ARRAY };

  // Gather large memory allocators in 1 place. This isn't the right
  // place, ultimately, but it is convenient! Having these available
  // allows the recursive sorts to share them.
  // AlignLadleSpace and PrefOfMem are in use only during a deal (and Close())
  // We must delay the recursion (which means the insertion sort)
  // until the semi sorted data has been returned to the source buffer
  // before BucketBuffer and MergedOfMem can be safely reused.
  AlignedBlock<true> BucketBuffer("BucketBuffer");
  AlignedBlock<true> AlignLadleSpace("AlignLadleSpace");
  // Overflow memory blocks.
  AlignedBlock<true> PrefOfMem("PrefOfMem");
  AlignedBlock<true> MergedOfMem("MergedOfMem");
  
// Support utilities directly connected with the sort infrastructure.

#include "DataInfo.H"
#include "DealByteAccess.H"
#include "PassInfoStruct.H"
#include "MuAligner.H"
  

#include "IlElem.H"
  template<typename SRec_T>
  using IlElemTemplate = IlElem<SRec_T, IlN>;
  
#include "GuessBkts.h"
#include "BitFieldAccess.H"

  typedef ChronoStopWatch<(TIMINGS > 1)> StopWatch;

  template<typename SRec_T>
  void DoSort(SRec_T *src_v, CountType src_n, int n_effective_flds);
  
  // Following includes expect ConstConfig context.

  //#include "LdlNdxManip.H"
#include "DealerBase.H"

  template<typename Elem_T>class InsertionSortPass;
  
  // Following select a version of each type of module.
#include "FlushMovers.H"
#include "Overflow.H"  // Not a module selector, but needs to fit here
#include "DealerCores.H"  
  //#include "LadleManager.H"  // Not a module selector, but needs to fit here
#include "DataPoolSet.h"
#include "FlusherSet.h"
#include "DealerExec.H"  // Not a module selector, but needs to fit here
#include "DataPoolPair.H"  // Not a module selector, but needs to fit here
#include "SdSet.h"

#include "InsertionSortPass.H"

  template<typename SRec_T>
  void DoSort(SRec_T *src_v, CountType src_n) {
    
    SdTemplate<SRec_T> SdObj/*(src_v, src_n)*/;

    SdObj.DoSort(src_v, src_n);
  }
  
};


#endif

