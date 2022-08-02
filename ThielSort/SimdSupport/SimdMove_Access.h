/* $Header:  $ */ 

/********************************************************************
 *
 *                            SimdMove_Access.H
 *
 * Access to components of Simd_Move_Type<>
 *
 * Since this may be included in more than one class, it doesn't use
 * the normal protection against multiple includes.
 *
 * Note that the "Simd_Move_Type" base class name follow my convention
 * for template parameter class names, so this list can be applied to
 * any variation on the SimdMove class.
 ********************************************************************/

using Simd_Move_Type::UtsPerSimd;
using Simd_Move_Type::UtsPerBlock;
using Simd_Move_Type::ModSimd;
using Simd_Move_Type::ModSimdBlk;
using Simd_Move_Type::RoundDownBlk;

using Simd_Move_Type::DoPrefetch;
using Simd_Move_Type::RevBulkMove;
using Simd_Move_Type::MoveSimds;
using Simd_Move_Type::MoveUts;
using Simd_Move_Type::MoveBlock;
using Simd_Move_Type::RevMoveBlock;


