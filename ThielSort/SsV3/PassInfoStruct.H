/* $Header:  $ */ 

#ifndef PassInfoStructH
#define PassInfoStructH

#include "ChronoStopWatch.H"

/********************************************************************
 *
 *                            PassInfoStruct.H
 *
 * Accumulate info about individual passes to be printed at the
 * end. Specifically includes functions to start and stop the
 * timer. If template parameter Timing_On is true the timing functions
 * are active. Otherwise, they are noops which should be ignored by an
 * optimizing compiler. Times are printed as microseconds.
 *
 * 2020-05-13 : Add support for intermittant timer. Defines functions
 * to turn timing on and off repeatedly to sum times in code that is
 * interspersed among code blocks we want to exclude. Original
 * functions give same results as before. Additional functions:
 *
 * Start() : Version without parameters to (re)start intermittant timer.
 *           Note that original Start(...) calls this one.
 * 
 * Stop() : Stops timer and adds time of current interval to
 *          Total. Synonym for End().
 *
 * 2020-08-01 : Add functions to return PassTimer and FlushTimer
 * times, so caller can sum and average them.
 *
 * 2021-3-16 : Control whether timing is active or not with a template
 *             parameter so we can allow total timings withoug pass
 *             timings. Check TIMINGS macro to support older code.
 ********************************************************************/

template<bool Timing_On = TIMINGS>class PassInfoStruct {
public:

  using StopWatch = ChronoStopWatch<Timing_On>;
  using ReportT = typename ChronoStopWatch<Timing_On>::ReportT;

  ChronoStopWatch<Timing_On> *GetFlushTimer() {
    //std::cout << "FlushTimer when addr taken:";
    //FlushTimer.Dump();
    return &FlushTimer;
  }
  ChronoStopWatch<Timing_On> *GetPassTimer() {
    return &PassTimer;
  }
  ReportT GetPassTotal() { return PassTimer.GetTotalTime(); }
  ReportT GetFlushTotal() { return FlushTimer.GetTotalTime(); }

  /**
   * Start this pass timer and set digit_number as its digit. 
   **/
  void Start(int digit_number) { Start("", digit_number); }

  /**
   * Init but don't start this timer.
   **/
  void Setup(char const *label_p) { DoSetup(label_p, -1); }
  
  /**
   * Start this function timer and set label_p as its title.
   **/
  void Start(char const *label_p) { Start(label_p, -1); }
  
  void Start() { PassTimer.Start(); }
  void End() { PassTimer.Stop(); }
  void Stop() { PassTimer.Stop(); }
  
  NO_INLINE void PrintPass(int pass) {
    if (!Timing_On) return;
    std::cout << "   Pass " << pass;
    std::cout << " digit " << DigitNumber;
    //std::cout << " " << DcLabel;  // Ignore DcLabel.
    ReportT passTime =  PassTimer.GetTotalTime();
    ReportT flushTime = FlushTimer.GetTotalTime();
    PassTimer.PrintLabeled(": deal=", passTime - flushTime);
    PassTimer.PrintLabeled(", flush=", flushTime);
    PassTimer.PrintLabeled(", total=", passTime);
    std::cout << '\n';
  }
  void PrintLabeled() {
    PassTimer.PrintTotalTime(DcLabel);
  }
  // Ignore DcLabel and use label_p.
  void PrintLabeled(char const *label_p) {
    PassTimer.PrintTotalTime(label_p);
  }
  void PrintLabeled(char const *label_p, ReportT time_val) {
    PassTimer.PrintLabeled(label_p, time_val);
  }
  // Debug. Dump values without changing anything
  void Dump(int ndx, const char *dump_label = nullptr) {
    std::cout << "Dump PassInfoStruct[" << ndx << "]:";
    if (dump_label) std::cout << "(label=" << dump_label << ") ";
    std::cout << "DcLabel: " << ((DcLabel) ? DcLabel : "null") <<
      ", DigitNumber=" <<
      DigitNumber << std::endl;
    std::cout << "PassTimer: "; PassTimer.Dump();
    std::cout << "FlushTimer: "; FlushTimer.Dump();
  }
private:

  ChronoStopWatch<Timing_On> PassTimer{};
  ChronoStopWatch<Timing_On> FlushTimer{};
  char const *DcLabel;
  int DigitNumber;

  void Start(char const *dc_label, int digit_number) {
    DoSetup(dc_label, digit_number);
    Start();
  }

  void DoSetup(char const *dc_label, int digit_number) {
    if (!Timing_On) return;
    DcLabel = dc_label;
    DigitNumber = digit_number;
  }
};

#endif
