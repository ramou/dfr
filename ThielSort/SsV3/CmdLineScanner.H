/* $Header:  $ */ 
#include <string>
#include "iostream"
#include <sstream>

#ifndef CmdLineScannerH
#define CmdLineScannerH

/********************************************************************
 *
 *                            CmdLineScanner.H
 *
 *
 ********************************************************************/

// CmdLineScanner                                           CmdLineScanner
class CmdLineScanner {
  const int argc;  
  const char **const argv;
  // Next item to look at in argv. argv[ArgNdx] is the "current arg".
  int ArgNdx;
public:
  // Default skips executable name and start with first param.
  CmdLineScanner(int arg_c, const char **arg_v, int arg_ndx = 1) :
    argc(arg_c), argv(arg_v), ArgNdx(arg_ndx) {}
    
  // If end of aragv (ArgNdx >= argc) return true. Else false.
  bool IsEnd() { return (ArgNdx >= argc); }
  int GetNdx() { return ArgNdx; }

  // NOTE: For all functions that return info from the arg list, we
  // move to the next arg when we find (and return) an expected value,
  // and don't move to the next arg when we don't find an expected
  // value.
  
  // Expect an unsigned int value at argv[ArgNdx]. If we find one,
  // store it into *target, advance ArgNdx, and return true. Else
  // return false.
  template<typename ElemT>
  bool ExpectUnsignedValue(ElemT *target) {
    if (ArgNdx >= argc) return false;
    if (!std::isdigit(argv[ArgNdx][0])) return false;
    std::string s(argv[ArgNdx]);
    *target = std::stod(argv[ArgNdx],nullptr);
    ArgNdx++;
    return true;
  }
  /**
   * Expect the current arg to be "-x", where x is a single
   * character. If so, return the value of the char following the
   * minus (it will be 0 if the minus stands alone) and move to the
   * next arg.  
   *
   * If we are at the end of the arg list, or the arg doesn't start
   * with a '-', return -1.
   *
   * If we have "-x" and x is not a single char, return -2. Use
   * ExpectMinusStr() to get the arg.
   **/
  int ExpectMinus() {
    if (ArgNdx >= argc) return -1;
    if (argv[ArgNdx][0] != '-') return -1;
    if (argv[ArgNdx][1] == 0) return argv[ArgNdx++][1];  // Just a minus
    if (argv[ArgNdx][2] == 0) return argv[ArgNdx++][1];  // expected result!
    return -2;
  }
  /**
   * Expect current arg to be "-x".  If the arg begins with '-' return
   * the entire arg (even if the arg is only a '-' or has only a
   * single character) and move to the next arg. Otherwise, return
   * nullptr and don't move to the next arg. This is designed to be
   * used after ExpectMinus() return -2.
   **/  
  const char* ExpectMinusStr() {
    if (ArgNdx >= argc) return nullptr;
    if (argv[ArgNdx][0] != '-') return argv[ArgNdx++];
    return nullptr;
  }
  /**
   * Expect the current arg to be a string not starting with '-'. If
   * so, store a pointer to it into target, advance ArgNdx, and return
   * true. If not, do nothing and return false.
   **/
  bool ExpectString(const char** target) {
    if (ArgNdx >= argc) return false;
    if (argv[ArgNdx][0] == '-') return false;
    *target = argv[ArgNdx++];
    return true;
  }
};  // CmdLineScanner

#endif
