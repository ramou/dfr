/**
 * Allows us to see what c++ makes of a collection of include files
 * without dealing with the cascade of irrelevant errors that can come
 * from fr.hpp.
 **/

#include <cassert>
#include <iostream>
#include <bits/stdc++.h>

#define isominch

#include "OddsAndEnds.h"

#include "Config.h"

typedef unsigned long long EType;

#ifndef N_BITS_TO_SORT
#define N_BITS_TO_SORT 0
#endif

struct HelpStruct {
  const char Tag;
  const bool ExtraInfo;  // If true, execute extra code to give info.
  const char *Text;
};
HelpStruct Help[] = {
  {'v', false, "Sort and verify, no config dump"},
  {'s', false, "sort only: no verify, no config dump"},
  {'r', true,  "random number generator specification: ..."},
  {'n', false, "option is n to sort"},
  {'h', false, "print this message and exit"},
  {'\0', true, ""},
};
void PrintExtra(char tag) {
  switch (tag) {
  case 'r': { RandomTool::CmdLineHelper(); break; }
  }
}
void PrintHelp() {
  std::cout << "Arg list is [arg] [arg] ... [nnnnnn]\n";
  std::cout << "arg is -? [options], where ? is a letter and options"
    " (if present) provides extra info.\n";
  std::cout << "If present, nnnnnn"
    " is the last option and is the number of values to sort.\n";
  for (int i = 0; Help[i].Tag; i++) {  // Loop til end encountered
    std::cout << '-' << Help[i].Tag << ": " << Help[i].Text << '\n';
    if (Help[i].ExtraInfo) PrintExtra(Help[i].Tag);
  }
}
int main(int argc, const char *argv[]) {
  AlignedBlock<> ValuesAB;
  bool Verify = true;
  CmdLineScanner cls(argc, argv);

  RandomTool::RandomBC<EType> *rng{ nullptr };

  long long n = 500027;  // Flag as default N.

  // Look for options: -? until first with no '-'
  // Default: DumpConfig, sort, verify.
  do {
    int tag = cls.ExpectMinus();
    if (tag < 0) break;  // End of args or not -...
    switch (tag) {
    case 'v': { ConstConfig::DumpConfig = false; Verify = true; break; }
    case 's': { ConstConfig::DumpConfig = false; Verify = false; break; }
    case 'r': { 
        rng = RandomTool:: template GetGenClass<EType>(&cls); 
	std::cout << "Random seed = " << rng->Seed() << std::endl;
	break;
    }
    case 'n': {
      if (!cls.ExpectUnsignedValue(&n)) ErrTag("Bad N", true);
      break;
    }
    case 'h': { PrintHelp(); exit(0); }
    default : {
      std::cout << "Unknown option is " << (char)tag << std::endl;
      exit(666);
    }
    }
  } while (true);
  if (rng == nullptr)  // Use default rng: Uniform, 64 bits unsigned.
    rng = RandomTool:: template GetGenClass<EType>(nullptr);
  
  rng->PrintConfig();

  cls.ExpectUnsignedValue(&n);  // n not changed unless number found.

  EType *Values; ValuesAB.GetAdjustedAdr(&Values, n);
  

  //auto seed = MakeRandomVec(Values, n, false, N_BITS_TO_SORT);
  //auto seed = MakeNormalRandomVec(Values, n, 1l<<32, 1<<30);
  rng->MakeRandomVec(Values, n);

  ConstConfig::DoSort(Values, n); 
  
#ifndef PROFILING
  if (Verify) {
    std::cout << "Starting verify" << std::endl;
    EType *DupVec = new EType[n];
    rng->MakeRandomVec(DupVec, n);
    //MakeNormalRandomVec(DupVec, n, 1l<<32, 1<<30, seed);
    //MakeSeededRandomVec(DupVec, n, seed, N_BITS_TO_SORT);
  
    std::sort(DupVec, DupVec + n);
    int errN = 0;
    for (int i = 0; i < n; i++) {
      if (DupVec[i] != Values[i]) {
	errN++;
	if (errN < 10) {
	  std::cout << "Error at index " << i << ", found "
		    << std::hex << Values[i] << " should be " << std::hex
		    << DupVec[i] << "\n" << std::dec;
	}
      }
    }
    if (errN > 0) {
      std::cout << "*****" << errN << " errors found.\n";
    } else {
      std::cout << "*****Sorted with no errors\n";
    }
  } else {
    std::cout << "*****Order not checked.\n";
  }
#endif
  
  return 0;
}
