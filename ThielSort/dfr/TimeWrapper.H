/* $Header:  $ */ 

#ifndef TimeWrapperH
#define TimeWrapperH

/********************************************************************
 *
 *                            TimeWrapper.H
 *
 * Wrap the PassInfoStruct vector, pass selector enum, and print
 * functions. Unclutters DataPoolPair and allows use for test cases
 * that don't want to use DataPoolPair (e.g., an Sd class that just
 * times std::sort).
 *
 * There are 3 kinds of timers (but all are stopwatch timers):
 *
 * - pass timers: They are allocated in order by calling
 *   StartPassTimer(), which starts them and sets them as the current
 *   pass timer. They are stopped by calling StopPassTimer().
 *
 * - function timers: They are allocated in order by calling
 *   StartFnTimer() and stopped by calling StopFnTimer(). When a
 *   function timer is started it becomes the current function timer,
 *   and StopFnTimer() with no parameter will stop it. 
 *
 * - default timers: TotalTimer and SetupTimer. TotalTimer is started
 *   by TimeWInit() and stopped by PrintInfo(). SetupTimer is also
 *   started by TimeWInit(), but must be stopped at the end of the
 *   setup phase by calling StopSetup().
 *
 * Notes:
 *
 * - an initial pass (e.g., to sample the live bits or digit counts)
 *   can be timed as a function or just billed to setup time. 
 *
 * - The max number of passes is MaxNPasses. The function timers
 *   (i.e., CopyTimeNdx) are in addition to that.
 *   
 * - The function timers are actually the pass timers of "extra"
 *   objects of PassTimerT (in vector PassInfo) to allow a uniform
 *   interface to all timers.
 *
 * - *** Call TimeWInit() before any other functions ***
 *
 * - There can be only 1 current function timer at any instant. A new
 *   call to StartFnTimer() before the current one is stopped doesn't
 *   stop the old one, but replaces it as the current timer. To nest
 *   timers (e.g., to separately time parts of a function), save the
 *   timer id returned by StartFnTimer() and pass it to StopFnTimer()
 *   to stop that timer. The setup timer uses this mechanism
 *   automatically, so there can be multiple timed segments in the
 *   "setup", (but nested functions will be included in the setup
 *   time).
 *
 * - Similarly, there is only 1 current pass timer at any instant, but
 *   they can be stopped and started explicitly by remembering the
 *   timer id of the pass timers.
 ********************************************************************/

#include "PassInfoStruct.H"

// TimeWrapper                                            TimeWrapper
class TimeWrapper {

public:

  using TotalTimerT = PassInfoStruct<(TIMINGS > 0)>;
  using PassTimerT =  PassInfoStruct<(TIMINGS > 1)>;
  using StopWatch = PassTimerT::StopWatch;
  
  static const int MaxNPasses = 20;  // Becomes index for total time timer.
  static const int MaxDfltTimers = 1;
  static const int MaxNFunctions = 5;  
  static const int MaxTtlTimers = MaxNPasses+MaxDfltTimers+ MaxNFunctions;
  
  //////////////////////////////////////////////////////////////////////
  //            Constructors / Initialization / Destructor

  //TimeWrapper(){}
  //~TimeWrapper(){}
  // Init timers. Return the timer id for the setup timer. n_to_sort
  // != 0 appears on the line with total time.
  int TimeWInit(long n_to_sort = 0) {
    NToSort = n_to_sort;
    TotalTimer.Start("Total time = ");
    return SetupTid = StartFnTimer("   Setup time = ");
  }

  //////////////////////////////////////////////////////////////////////
  //                   Other functions.
  
  /**
   * Make next timer the current pass timer. Start it and set
   * digit_number > 0 as its digit for printing. Return the timer id
   * to allow nested and intermittent timers.
   **/
  int StartPassTimer(int digit_number) {
    PassId++;
    if (PassId >= MaxNPasses) MyAbort("Too many passes for TimeWrapper");
    //PassInfo[PassId].Dump(PassId, "Before initial start");
    PassInfo[PassId].Start(digit_number);
    //PassInfo[PassId].Dump(PassId, "After Init");
    return PassId;
  }
  StopWatch *GetFlushTimer(int tmr_id) {  // Current pass
    return PassInfo[tmr_id].GetFlushTimer();
  }
  void StopPassTimer() { PassInfo[PassId].Stop(); }
  /**
   * Allocate a new function timer, set its label to label_p, and
   * return the timer id. The timer isn't started.
   **/
  int AllocateFnTimer(const char *label_p) {
    if (label_p == nullptr) ErrTag("WTF");
    FunctionId++;
    if (FunctionId >= MaxTtlTimers)
      MyAbort("Too many functions for TimeWrapper");
    PassInfo[FunctionId].Setup(label_p);
    return FunctionId;
  }
  /**
   * Allocate a new function timer, set its label to label_p, and
   * return the timer id. If start_tmr is true, the timer is started.
   **/
  int StartFnTimer(const char *label_p, bool start_tmr = true) {
    int tmrId = AllocateFnTimer(label_p);
    if (start_tmr) PassInfo[tmrId].Start(label_p);
    return tmrId;
  }
  // Stop current function timer.
  void StopFnTimer() { PassInfo[FunctionId].Stop(); }
  // Stop a function timer. Deprecated: Use StopTimer()
  void StopFnTimer(int tmr_id) { StopTimer(tmr_id); }

  // Stop and restart timer tmr_id (any timer).
  void StopTimer(int tmr_id) { PassInfo[tmr_id].Stop(); }
  // (Re)start a possibly new fn timer.
  int ReStartTimer(int tmr_id, const char *label_p) {
    if (tmr_id < 0) tmr_id = StartFnTimer(label_p);
    else PassInfo[tmr_id].Start();
    return tmr_id;
  }
  // (Re)start a possibly new pass timer.
  int ReStartTimer(int tmr_id, int digit_number) {
    if (tmr_id < 0) tmr_id = StartPassTimer(digit_number);
    else PassInfo[tmr_id].Start();
    return tmr_id;
  }
  // Restart a previously used timer (either type).
  void ReStartTimer(int tmr_id) { PassInfo[tmr_id].Start(); }
  void StopSetup() { StopTimer(SetupTid); }
  
// PrintInfo                                                      PrintInfo  
  NO_INLINE void PrintInfo() {
    TotalTimer.End();
    if (TIMINGS > 1) {
      // Print pass timers
      for (int i = 0; i <= PassId; i++) {
	PassInfo[i].PrintPass(i);
      }
      /**
       * Pass 1 and the last pass times are often anomalous due to
       * initial use of new memory (Pass 1 times) and storage back to
       * original matrix if last pass stores directly. We compute
       * average values for passes 2 : PassId-1.
       **/
      if (PassId > 1) {
	using ReportT = PassTimerT::ReportT;
	ReportT flushTim = PassInfo[1].GetFlushTotal();
	ReportT passTim = PassInfo[1].GetPassTotal();
	for(int i = 2; i < PassId; i++) {
	  flushTim += PassInfo[i].GetFlushTotal();
	  passTim += PassInfo[i].GetPassTotal();
	}
	int np = PassId - 1;
	StopWatch::PrintLabeled("         Inner averages : deal= ",
				(passTim - flushTim) / np);
	StopWatch::PrintLabeled(", flush= ", flushTim / np);
	PutLn();
      }
      // print function timers
      for (int i = MaxNPasses; i <= FunctionId; i++) {
	PassInfo[i].PrintLabeled();
      }
    }
    std::cout << "\n";
    if (NToSort) {
      std::cout << "N = ";
      std::cout.width(10); std::cout << NToSort;
    }
    TotalTimer.PrintLabeled("   Total time= ");
    PutLn();
  }

private:
  long NToSort = 0;
  PassTimerT PassInfo[MaxTtlTimers];
  TotalTimerT TotalTimer;
  
  // Automatic counter for true pass timers. Identifies current pass.
  // Incremented on each call to StartPassTimer().
  // Incremented before use, thus init to -1 to start at PassInfo[0].
  int PassId = -1;
  // Automatic counter for function timers. Identifies current function.
  // Incremented on each call to StartFnTimer(). Includes default timers
  int FunctionId = MaxNPasses - 1;
  int SetupTid = -1;
};  // TimeWrapper

#endif
