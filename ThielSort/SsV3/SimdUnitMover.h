/* $Header:  $ */ 

#ifndef SimdUnitMoverH
#define SimdUnitMoverH

/********************************************************************
 *
 *                            SimdUnitMover.H
 *
 * The SimdUnitMover namespace contains support for using the simd
 * instructions to move multiple records at a time, using NT
 * instructions when indicated. The simd unit (SDU) is the 16, 32, or
 * 64 (not yet) byte vector selected using the SIMD_MOVE_SIZE
 * parameter defined in the makefile to reflect the current hardware.
 *
 * We define 2 versions of SDU that a program may use:
 *
 * LocalSdu : register variables and variables expected to be accessed
 *            via the hardware cache.
 *
 * RemoteSdu : memory variables (normally vectors) which should be
 *             accessed by non temporal (NT) instructions which bypass
 *             the hardware cache when possible.
 *
 * AVX limits using NT instructions to storing data only.
 *
 * Assignment (=) operators are provided to use the correct
 * instructions when doing assignments within and between these types.
 *
 * At this time only block moves of these types is supported here.
 ********************************************************************/

#include "SrtSimdBase.H"

namespace SimdUnitMover {

  /**
   * Constants for templates to specify LocalSdu (use normal
   * instructions to access) or RemoteSdu (use NT instructions when
   * possible).
   **/
  enum MemoryType { LocalMem, RemoteMem };
  
  /**
   * Template for the basic SDU. 
   **/  
  template<MemoryType Mem_Type, int Move_Size = SIMD_MOVE_SIZE>
    struct SimdUnit : public SrtSimdBase<Move_Size> {

    typedef SrtSimdBase<Move_Size> BaseT;
    typedef typename BaseT::IntelSdType VecT;
    static const int UnitLen = 1 << Move_Size;

    SimdUnit &operator= (const BaseT &v) {
      if (Mem_Type == LocalMem) this->MoveIn(v);
      else this->NtMoveIn(v);
      return *this;
    }
    SimdUnit &operator= (SimdUnit<LocalMem, Move_Size> &v) {
      if (Mem_Type == LocalMem) this->MoveIn((SrtSimdBase<Move_Size>)v);
      else this->NtMoveIn((SrtSimdBase<Move_Size>)v);
      return *this;
    }
    SimdUnit &operator= (SimdUnit<RemoteMem, Move_Size> &v) {
      if (Mem_Type == LocalMem) this->MoveIn((SrtSimdBase<Move_Size>)v);
      else this->NtMoveIn((SrtSimdBase<Move_Size>)v);
      return *this;
    }
    
    // Convert SortRecType* to SimdUnit*. Function makes code less messy.
    static SimdUnit *SimdPtr(void *p)
    { return reinterpret_cast<SimdUnit *>(p); }
    static const SimdUnit *SimdPtr(const void *p)
    { return reinterpret_cast<const SimdUnit *>(p); }
    
  }__attribute__ ((aligned(UnitLen)));  // SimdUnit< ? >

  template<int Move_Size> using LocalSdu = SimdUnit<LocalMem, Move_Size>;
  template<int Move_Size> using RemoteSdu = SimdUnit<RemoteMem, Move_Size>;
  
};  // SimdUnitMover

#ifdef SimdUnitMover_home
static const RegisterSource SimdUnitMover_H_Registry("$Header:  $");
#endif

#endif
