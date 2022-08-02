/* $Header:  $ */ 

#ifndef VectorTypeH
#define VectorTypeH

/********************************************************************
 *
 *                            VectorType.H
 *
 * A templated class to serve as the base type for any simd vector
 * declared by the user. This file includes a bunch of other files to
 * provide the nuts and bolts, but class VectorType<> provides the
 * complete user interface and is the only file that needs to be
 * included directly. It implements a generic vector type that is
 * specialized by template arguments that a select simd unit size and
 * an element type from the lists below.
 *
 ********************************************************************/

#include "IntrinSetup.h"
#include "SimdDcls.h"
#include "VTBase.H"
#include "TypeSpecific.h"

namespace SimdSupport {

// VectorType                                                     VectorType
  template<typename Elem_T, int Log_Su_Len>
  class VectorType : public VTBase<Log_Su_Len>,
		     private ElemBase<Elem_T, Log_Su_Len> {
    typedef VTBase<Log_Su_Len> VtBc;  // VectorType base class.
    // Imp for ops, access to intrinsics:
    typedef ElemBase<Elem_T, Log_Su_Len> VtImp; 

  public:
    typedef typename VtImp::VecType EVecType;  // vec of Elem_T

    using VtBc::SuLen;
    // How many Elem_T s are in an EVecType?
    static const int NPerSimd = SuLen / sizeof(Elem_T);
    
    //////////////////////////////////////////////////////////////////////
    //            Constructors / Initialization / Destructor

    using VtBc::VtBc;
    VectorType(Elem_T v) { *this = BroadCast(v); }
    VectorType() :VtBc() {}
    VectorType(EVecType &src) { VtBc::Store(src); }

    //~VectorType();

    //////////////////////////////////////////////////////////////////////
    //                   Other functions.
    using VtImp::BroadCast;
    using VtImp::Extract;
    using VtImp::Insert;
    using VtImp::ShiftLeft;
    using VtImp::ShiftRight;
    using VtImp::IVMult;
    using VtImp::IVAdd;
    using VtImp::IVSub;
    using VtImp::AndInstr;
    using VtImp::AndNotInstr;
    using VtImp::OrInstr;
    using VtImp::XorInstr;
      
    VectorType &operator =(const VtBc &src)
    { VtBc::Store(src); return *this; }
    VectorType &operator =(const EVecType &src)
    { VtBc::Store(src); return *this; }
    // Return element in slot ndx. Not a ref to store to.
    Elem_T operator[](int ndx) { return GetElem(*this, ndx); }
    
    void DumpVec(const char *title) const { VtImp::DumpVec(*this, title); }
    /**
     * Reset the vector to a new value:
     *
     * if n >= NPerSimd, set the vector to the first NPerSimd elements
     * of init_v.  If n < NPerSimd, set the first n elements of the
     * vector to the corr elements of init_v, and the rest to 0.
     **/    
    void Load(Elem_T *init_v, int n) {
      if (n >= NPerSimd) n = NPerSimd;
      else *this = XorInstr(*this, *this);  // ~[n:NPerSimd-1] = 0 
      for (int i = 0; i < n; i++) *this = Insert(*this, i, init_v[i]);
    }
    // Arithmetic
    VectorType operator +(const VtBc &right_op) {
      return IVAdd(*this, right_op);
    }
    VectorType operator -(const VtBc &right_op) {
      return IVSub(*this, right_op);
    }
    // Casting
    static int IsAligned(const void *p) {
      return ((unsigned long)p & (sizeof(VtBc) - 1)) == 0;
    }
    static VectorType *CastAddr(void *p)
    { return reinterpret_cast<VectorType*>(p); }
    static const VectorType *CastAddr(const void *p)
    { return reinterpret_cast<const VectorType*>(p); }
    // Shift and logical
    VectorType operator <<(int bit_cnt) { return ShiftLeft(*this, bit_cnt); }
    VectorType operator >>(int bit_cnt) { return ShiftRight(*this, bit_cnt); }
    VectorType operator &(const VtBc &right_op)
    { return AndInstr(*this, right_op); }
    VectorType operator |(const VtBc &right_op)
    { return OrInstr(*this, right_op); }
    VectorType operator ^(const VtBc &right_op)
    { return XorInstr(*this, right_op); }
    VectorType AndNot(const VtBc &right_op)  // *this &~ right_op
    { return AndNotInstr(right_op, *this); }
  private:
    Elem_T GetElem(const EVecType &vec, int ndx) {
      return vec[ndx];
    }
  } SIMD_ALIGN(1<<Log_Su_Len);  // VectorType

};  // namespace SimdSupport
#ifdef VectorType_home
static const RegisterSource VectorType_H_Registry("$Header:  $");
#endif

#endif
