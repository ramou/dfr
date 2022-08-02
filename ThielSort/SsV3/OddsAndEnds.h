/* $Header:  $ */ 

#ifndef OddsAndEndsH
#define OddsAndEndsH

#include <cstdint>
#include "iomanip"

/********************************************************************
 *
 *                            OddsAndEnds.H
 *
 ********************************************************************/

#define NO_INLINE __attribute__ ((__noinline__))
#define INLINE_ATT __attribute__ ((__always_inline__))

typedef unsigned char UByte;  // 0:255
typedef signed char SByte;  // -127:127
typedef unsigned short UShort;
typedef uint32_t UInt;
typedef uint64_t ULong;  
#if 0
typedef unsigned __int64 uint64_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int8 uint8_t;
#endif

// Give us cheap access to function / block on first call. Usually to
// disassemble the code. No calling code generated unless GDBHOOK defined.
NO_INLINE void RealGdbHook(int caller_id) {
  std::cout << "GdbHook(" << caller_id << ")\n";
}
/**
 * Print to signal if this gets left in by mistake, but don't print
 * too often!
 **/
INLINE_ATT void GdbHook(int caller_id) {
#ifdef GDBHOOK
  static long Flags = 0;
  int flg = 1 << caller_id;
  if ((Flags & flg) == 0) {  // First call
    Flags |= flg;
    RealGdbHook(caller_id);
  }
#endif
}

// Return pwr = the number of low order binary 0 in representation of
// i. Assuming i is a power of 2, (1 << pwr) == i.
constexpr int FindPower(int i) {
  if (i & 1) return 0;
  else return 1 + FindPower(i >> 1);
}

void ErrTag(const char *msg = nullptr, bool die = false) {
  if (msg == nullptr) msg = "No Message";
  std::cout << msg << "\n";
  if (die) exit(91);
}
[[nreturn]] void MyAbort(const char *msg) {
  std::cout << msg << std::endl;
  exit(999);
}

// Helpers for ConfirmTargetBkts()
void InBrackets(int b) { std::cout << '[' << b << ']'; }
void BadDealByte(int b, const char *p_o_ptr, int ndx, int val) {
  std::cout << p_o_ptr;
  InBrackets(b);
  std::cout <<  " item";
  InBrackets(ndx);
  std::cout <<  " has bad dealbyte: " << val << std::endl;
}

// Helper for DumpIntVec()
void DumpIntLine(const int *vec, int start_ndx, int n, int fld_width) {
  std::cout << std::setw(2) << start_ndx << ": ";
  vec += start_ndx;
  n--;
  for (int i = 0; i <= n; i++)
    std::cout << std::setw(fld_width) << vec[i] <<
      ((i < n)?", " : "\n");
}
/**
 * Dump vec[0:n-1] to stdout. Format:
 * label_p
 * x: vec[x], vec[x+1], ...
 * ...
 *
 * x is 2 or more digits, vec[?] is fld_width or more digits. Number
 * of items per line is LineLen / (fld_width+2).
 **/
void DumpIntVec(const int *vec, int n, const char *label_p,
		int fld_width = 4) {
  static const int LineLen = 60;
  int nPerLine = LineLen / (fld_width+2);
  std::cout << label_p << '\n';
  int end = n - (nPerLine - 1);
  int lineStart = 0;  // index in vec of first in next print line.
  for (; lineStart <=end; lineStart += nPerLine) {
    DumpIntLine(vec, lineStart, nPerLine, fld_width);
  }
  if (lineStart < end) {  // need short line to finish
    DumpIntLine(vec, lineStart, nPerLine - lineStart, fld_width);
  }
}

// IO helpers
void PutLn() { std::cout << '\n'; }
void PutInt(long val, int width) {
  std::cout.width(width); std::cout << val;
}

/**
 * This simple sort checks N and does an insertion sort if it is
 * smaller than IsLimit, and calls std::sort() otherwise.
 * If n < 2, just return.
 **/
template<typename T>void MiniSort(T* src_v, int n) {
  constexpr int IsLimit = 15;
  if (n >= IsLimit) {
    std::sort(src_v, src_v + n);
    return;
  }
  if (n < 2) return;
  T smallest = src_v[0];  // Put in register if T is a simple number.
  for (int i = 1; i < n; i++) {
    T curVal = src_v[i];
    if (curVal < smallest) {
      for (int j = i; j > 0; j--) src_v[j] = src_v[j-1];
      src_v[0] = smallest = curVal;
    } else {  // Sift curVal down
      T *p = src_v + i - 1;
      while (p[0] > curVal) { p[1] = p[0]; p--; }
      p[1] = curVal;
    }
  }
}

#ifndef CACHE_LINE_LEN
#define CACHE_LINE_LEN 64
#endif
static const int CacheLineLen = CACHE_LINE_LEN;


#endif
