// -*- C++ -*-
/***************************************************************************
 * blitz/meta/vecassign.h   TinyVector assignment metaprogram
 *
 * $Id$
 *
 * Copyright (C) 1997-2011 Todd Veldhuizen <tveldhui@acm.org>
 *
 * This file is a part of Blitz.
 *
 * Blitz is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * Blitz is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public 
 * License along with Blitz.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Suggestions:          blitz-devel@lists.sourceforge.net
 * Bugs:                 blitz-support@lists.sourceforge.net    
 *
 * For more information, please see the Blitz++ Home Page:
 *    https://sourceforge.net/projects/blitz/
 *
 ***************************************************************************/

#ifndef BZ_META_VECASSIGN_H
#define BZ_META_VECASSIGN_H

BZ_NAMESPACE(blitz)

template<int N, int I, typename T_numtype> 
class _bz_meta_vecAssign {
public:
    // Calculate the SIMD vector width in T_numtype units.
    static const int vecwidth=
      BZ_SIMD_WIDTH>sizeof(T_numtype) ? BZ_SIMD_WIDTH/sizeof(T_numtype) : 1;
    static const int blocksize = N-I>=vecwidth ? vecwidth : 1;
    static const int fastLoopFlag = (I < N-blocksize) ? 1 : 0;
    static const int loopFlag = (I < N-1) ? 1 : 0;

    template<typename T_expr, typename T_updater>
    static inline void fastAssign(T_numtype* restrict data, T_expr expr, 
				  T_updater u)
    {
      // This metaprogram used to unroll the loops completely, but
      // this inhibits vector optimizations on icpc. If a SIMD width
      // is set, we therefore make a loop of the appropriate width to
      // facilitate vectorization. If BZ_SIMD_WIDTH==1, which is the
      // default, it will still unroll completely.

      // These pragmas are for icpc. Ideally this should be
      // compiler-independent, but currently icpc is the only compiler
      // for which they are used.
      #pragma ivdep
      #pragma vector aligned
      for(int i=I;i<I+blocksize;++i)
	u.update(data[i], expr._bz_fastAccess(i));

      _bz_meta_vecAssign<N * fastLoopFlag,
			 (I+blocksize) * fastLoopFlag, 
			 T_numtype>::fastAssign(data,expr,u);
    }

    template<typename T_vector, typename T_expr, typename T_updater>
    static inline void assign(T_vector& vec, T_expr expr, T_updater u)
    {
        u.update(vec[I], expr[I]);
        _bz_meta_vecAssign<N * loopFlag, (I+1) * loopFlag, T_numtype>
           ::assign(vec,expr,u);
    }

    template<typename T_vector, typename T_numtype2, typename T_updater>
    static inline void assignWithArgs(T_vector& vec, T_updater u,
        T_numtype2 x0, T_numtype2 x1=0, T_numtype2 x2=0, T_numtype2 x3=0,
        T_numtype2 x4=0, T_numtype2 x5=0, T_numtype2 x6=0, T_numtype2 x7=0,
        T_numtype2 x8=0, T_numtype2 x9=0)
    {
        u.update(vec[I], x0);
        _bz_meta_vecAssign<N * loopFlag, (I+1) * loopFlag, T_numtype>
            ::assignWithArgs(vec, u, x1, x2, x3, x4, x5, x6, x7, x8, x9);
    }
        
};

template<typename T_numtype>
class _bz_meta_vecAssign<0,0, T_numtype> {
public:
    template<typename T_expr, typename T_updater>
    static inline void fastAssign(T_numtype* restrict, T_expr, T_updater)
    { }

    template<typename T_vector, typename T_expr, typename T_updater>
    static inline void assign(T_vector&, T_expr, T_updater)
    { }

    template<typename T_vector, typename T_numtype2, typename T_updater>
    static inline void assignWithArgs(T_vector&, T_updater,
        T_numtype2, T_numtype2 =0, T_numtype2 =0, T_numtype2 =0,
        T_numtype2 =0, T_numtype2 =0, T_numtype2 =0, T_numtype2 =0,
        T_numtype2 =0, T_numtype2 =0)
    {
    }
};

BZ_NAMESPACE_END

#endif // BZ_META_ASSIGN_H
