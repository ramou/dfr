/* $Header:  $ */ 

#ifndef DataPoolSetH
#define DataPoolSetH

/********************************************************************
 *
 *                            DataPoolSet.h
 *
 * Collect the available DataPool .H files here and select 1 according
 * to the SELECT_DATAPOOL define constant.
 ********************************************************************/

// Generic datapool class template. Doesn't exist.
template<typename SRec_T, int Selector_Id>class DataPoolClass;

// Various specializations of DataPoolClass:

#include "DataPool.H"
//#include "LfDataPool.H"
//#include "LfV2DataPool.H"

template<typename SRec_T>
using DataPoolTemplate = DataPoolClass<SRec_T, SELECT_DATAPOOL>;

#endif
