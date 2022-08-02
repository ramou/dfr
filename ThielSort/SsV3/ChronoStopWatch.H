/* $Header:  $ */ 

#ifndef ChronoStopWatchH
#define ChronoStopWatchH

#include <cmath>
#include <chrono>

/********************************************************************
 *
 *                            ChronoStopWatch.H
 *
 * Measure elapsed time. Will accumulate multiple intervals. Format
 * the time into an internal buffer from which it can be printed.
 * Uses std::chrono::high_resolution_clock for timing. 
 * 
 * Times can be obtained for the sum of all intervals completed to
 * date (total time) or for the currently running (or most recently
 * stopped) interval (interval time). 
 *
 * By default, printing is done in microseconds. If template parameter
 * Timing is false, all action will be skipped, and an optimizing
 * compiler should make all calls EXCEPT gets and digit selectors,
 * which return values, disappear.
 *
 * If Count_Intervals is true, the default printer will include the
 * number of intervals that contributed to the time if it is > 1.
 *
 * NOTE: Not all options tested!
 ********************************************************************/

// CswBaseClass                                             CswBaseClass
/**
 * Untemplated base class with virtual functions to give us a handle
 * to create a ClockVector to simplify printing a set of timers. Only
 * contains private (to ChronoStopWatch) data members and functions.
 **/
class CswBaseClass {
public:
  using duration = std::chrono::high_resolution_clock::duration;
  using time_point = std::chrono::high_resolution_clock::time_point;
  
  // Support printing
  static const int MicroSecWidth = 6;   
  static const int NanoSecWidth = 9;   
  bool TimingIsOn;  // So we can see it outside.

  virtual
  void PrintTotalTime(char const *label_p, int n_frac_digits = MicroSecWidth) {
    std::cout << "CswBaseClass::PrintTotalTime called???\n";
    exit(0);
  }

  virtual
  void PrintIntervalTime
  (char const *label_p, int n_frac_digits = MicroSecWidth) {
    std::cout << "CswBaseClass::PrintIntervalTime called???\n";
    exit(0);
  }
protected:
  time_point StartTime{};
  duration Total{};
  int NIntervals = 0;
  
  time_point now() { return std::chrono::high_resolution_clock::now(); }

  CswBaseClass(int timing_is_on = true) : TimingIsOn(timing_is_on) {}
  
};
// ChronoStopWatch                                         ChronoStopWatch
template<bool Timing = true, bool Count_Intervals = false>
class ChronoStopWatch : public CswBaseClass {
public:
  using CswBaseClass::duration;
  using CswBaseClass::time_point;
  
  typedef long long int ReportT;  // type for reporting duration counts
  using DefltUnitT = std::chrono::duration<ReportT, std::micro>;

  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  // Constructor starts clock iff start_running == true (default false).
  ChronoStopWatch(bool start_running = false) : CswBaseClass(Timing) {
    if (start_running) Start();
  }
  //~ChronoStopWatch();

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.

  void Start() {
    if (Timing) {
      StartTime = now();
      if (Count_Intervals) NIntervals++;
    }
  }
  void Stop() {
    if (Timing) {
      //time_point startCopy{StartTime};
      Total += now() - StartTime;
      //if (Total.count() < 0) GdbHook(43);
      //std::cout << "Before Stop(): "; Dump();
      //std::cout << "After Stop(): "; Dump();
    }
  }
  
  /**
   * Reset the stopwatch. All times to date will be lost.
   **/
  void Reset() {
    if (Timing) {
      Total = duration{}; StartTime = time_point{}; NIntervals = 0;
    }
  }
  /**
   * Return the duration (total or interval time, as requested) as a
   * count of intervals of unit_size given by a std:: ratio name
   * (e.g., std::micro). Default is micro, so return value in that
   * case is microseconds. A currently running interval is not
   * included in GetTotalTime().
   **/
  template<typename UnitT>
  ReportT GetTotalTime(UnitT unit_size) {
    return std::chrono::duration_cast<UnitT>(Total).count();
  }
  ReportT GetTotalTime() {
    using UnitT = std::chrono::duration
      <ReportT, std::micro>;
    return std::chrono::duration_cast<UnitT>(Total).count();
  }
  template<typename UnitT>
  ReportT GetIntervalTime(UnitT unit_size) {
    return std::chrono::duration_cast<UnitT>(now() - StartTime).count();
  }
  ReportT GetIntervalTime() {
    using UnitT = std::chrono::duration
      <ReportT, std::micro>;
    return std::chrono::duration_cast<UnitT>(now() - StartTime).count();
  }
  /**
   * Return low order or high order digits of a ReportT value by using
   * fraction_part with mod or divide. Minor helpers for formatting
   * output.
   **/
  static ReportT LowDigits(ReportT val, ReportT fraction_part) {
    return val % fraction_part;
  }
  static ReportT HighDigits(ReportT val, ReportT fraction_part) {
    return val / fraction_part;
  }
  /**
   * Print "label_p space time_val[ in n intervals]",
   *        where " space " means 1 space
   *        and time_val is n...n if time_val < 10**n_frac_digits, 
   *                     or s.fraction otherwise, where fraction has
   *                        n_frac_digits digits.
   *        If the n_intervals parameter to PrintLabeled() > 1, the
   *        optional phrase (in [] above) will show it.
   *        No endl.
   **/
  void PrintTotalTime(char const *label_p, int n_frac_digits = MicroSecWidth) {
    PrintLabeled(label_p, GetTotalTime(), n_frac_digits,
		 Count_Intervals ? NIntervals : 1);
  }
  void PrintIntervalTime(char const *label_p,
			 int n_frac_digits = MicroSecWidth) {
    PrintLabeled(label_p, GetIntervalTime(), n_frac_digits);
  }  
  static void PrintLabeled(char const *label_p, ReportT time_val,
			   int n_frac_digits = MicroSecWidth,
			   int n_intervals = 1) {
    if (!Timing) return;
    if (label_p == nullptr) label_p = "";
    std::cout << label_p << ' ';
    ReportT unitsPerSec = pow(10, n_frac_digits);
    if (time_val > unitsPerSec) {  // print as secs.fraction
      std::cout << HighDigits(time_val, unitsPerSec) << '.';
      std::cout.width(n_frac_digits);
      char prev = std::cout.fill('0');
      std::cout << LowDigits(time_val, unitsPerSec);
      std::cout.fill(prev);
    } else {  // print as an int
      std::cout << time_val;
    }
    if (n_intervals > 1) {
      std::cout << "  in " << n_intervals << " intervals";
    }
  }
  /**
   * Since the print functions don't print eoln, the calling function
   * must. But the calling function may not know if timing (and thus
   * printing) is suppressed. PrintLn() will print an eoln iff timing
   * (and printing) is active.
   **/
  static void PrintLn() { if (Timing) std::cout << '\n'; }
  // Debug: dump local values without changing anything.
  void Dump() {
    time_point tt{};
    std::cout
      << "StartTime=" << std::chrono::duration_cast<DefltUnitT>
      (StartTime - time_point{}).count() <<
      ", Current total=" << GetTotalTime() << std::endl;
  }
};  // ChronoStopWatch

#if 0
#include <iostream>
using std::cout;

constexpr int N = 100000000;
void doubleme(int *p) {
  for (int i = 0; i < N; i++)  p[i] = p[i] + p[i];
}
int values[N];
int main() {
  ChronoStopWatch<true> t0;
  ChronoStopWatch<true> t1;
  ChronoStopWatch<true> t2;

  for (int i = 0; i < N; i++)  values[i] = i;
  
  t0.Start();
  t1.Start();
  doubleme(values);
  t1.Stop();
  t2.Start();
  doubleme(values);
  t2.Stop();
  t0.Stop();
  cout << "t1 = " << t1.GetTotalTime() << " t2= " << t2.GetTotalTime() << '\n';
  cout << "sum =" << t1.GetTotalTime() + t2.GetTotalTime() << '\n';
  cout << "t0 = " << t0.GetTotalTime() << '\n';


  
  return 0;
}
#endif
#endif
