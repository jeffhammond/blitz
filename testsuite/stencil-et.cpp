#include "testsuite.h"
#include <blitz/array.h>
#include <blitz/array/stencil-et.h>
#include <blitz/tinyvec-et.h>
#include <blitz/matrix.h>
#include <blitz/tinymatexpr.h>
#include <random/uniform.h>

BZ_USING_NAMESPACE(blitz)

// Tests that the various stencil operators work with
// expressions. Does NOT test that the stencils produce the correct
// output, only that they compile and that the output is consistent
// between different identical expressions.

typedef blitz::Array<double,1> array_1;
typedef blitz::Array<double,2> array_2;
typedef blitz::Array<double,3> array_3;
typedef Array<TinyVector<double, 2>, 2> array_2v;
typedef Array<TinyMatrix<double, 2, 2>, 2> array_2m;
typedef Array<TinyVector<double, 3>, 3> array_3v;
typedef Array<TinyMatrix<double, 3, 3>, 3> array_3m;

// test with functors
class doubler {
public:
  double operator()(double x) {return 2.0*x;}
  BZ_DECLARE_FUNCTOR(doubler);
};

class multiplier {
public:
  double operator()(double a, double b) {return a*b;}
  BZ_DECLARE_FUNCTOR2(multiplier);
};

/*
// Test two expressions for equality
template<typename T1, typename T2>
void test_expr(const T1& d1, const T2& d2)
{
  BZTEST(all(d1==d2));
}
*/
#define test_expr(d1,d2) BZTEST(all((d1)==(d2)));

// Test two vector expressions for equality
template<typename T1, typename T2>
void test_vexpr(const T1& d1, const T2& d2)
{
  Array<typename T1::T_numtype, T1::rank> a(d1),b(d2);
  for(int i=0; i<T1::T_numtype::numElements; ++i)
    BZTEST(all(a[i]==b[i]));
}

// Test two matrix expressions for equality
template<typename T1, typename T2>
void test_mexpr(const T1& d1, const T2& d2)
{
  // there appears to be no way to tell the size of a TinyMatrix... or to compare them.
  // so for now we're happy just to have the call succeed.
}


int main()
{
  // create some arrays to operate on
  const int sz=5;

  ranlib::Uniform<float> rnd;
  rnd.seed(42);

  array_3 field3(sz,sz,sz), result3(sz,sz,sz),
    fx(field3.shape()), fy(field3.shape()), fz(field3.shape());
  for(int i=0; i<field3.size();++i) {
    field3.data()[i]=rnd.random();
    fx.data()[i]=rnd.random();
    fy.data()[i]=rnd.random();
    fz.data()[i]=rnd.random();
  }
  array_2 field2(sz,sz), result2(field2.shape());
  field2=sin(0.5*(tensor::i+tensor::j));
  array_3v vfield3(field3.shape());
  vfield3[0]=fx;
  vfield3[1]=fy;
  vfield3[2]=fz;
  array_2v vfield2(sz,sz);
  vfield2[0]=vfield3(Range::all(), Range::all(), 0)[0];
  vfield2[1]=vfield3(Range::all(), Range::all(), 0)[1];

  doubler doubleit;
  multiplier multiplyit;

  // Now apply "all" possible stencil types to arrays and expressions,
  // as well as recursive applications

  // defined with BZ_ET_STENCIL:
  test_expr(Laplacian2D(field2), Laplacian2D(1.0*field2));
  test_expr(Laplacian2D(field2), 1.0*Laplacian2D(field2));
  test_expr(Laplacian2D(const_cast<const array_2&>(field2)),
	    Laplacian2D(1.0*field2));
  test_expr(Laplacian2D(Laplacian2D(field2)), 
	    Laplacian2D(Laplacian2D(1.0*field2)));
  test_expr(Laplacian2D(Laplacian2D(field2)), 
	    Laplacian2D(1.0*Laplacian2D(field2)));
  test_expr(Laplacian2D(Laplacian2D(field2+field2)), 
	    Laplacian2D(Laplacian2D(field2)+Laplacian2D(field2)));
  
  // and some more complicated expressions and assignments
  result2(_bz_shrinkDomain(result2.domain(),shape(-1,-1),shape(1,1))) =
    Laplacian2D(where(field2>0.5,0.,1.));
  test_expr(result2(_bz_shrinkDomain(result2.domain(),shape(-1,-1),shape(1,1))), 
	    Laplacian2D(where(field2>0.5,0.,1.)));
  test_expr(where(Laplacian2D(field2)>0.5,0.,1.),
	    where(Laplacian2D(2*field2)>1,0.,1.));
  test_expr(where(Laplacian2D(field2)>0.5, 
		  0.0*field2(_bz_shrinkDomain(result2.domain(),shape(-1,-1),shape(1,1))), 
		  0.0*field2(_bz_shrinkDomain(result2.domain(),shape(-1,-1),shape(1,1)))+2.0),
	    2*where(Laplacian2D(2*field2)>1.0, 0., 1.));
  test_expr(Laplacian2D(2.0*field2), Laplacian2D(doubleit(field2)));
  test_expr(Laplacian2D(field2*field3(0,Range::all(), Range::all())), 
	    Laplacian2D(multiplyit(field2, field3(0, Range::all(), Range::all()))));

  // defined with BZ_ET_STENCIL2:
  test_expr(div(vfield2[0],vfield2[1]),
	    div(vfield2[0],1.0*vfield2[1]));
  test_expr(div(vfield2[0],vfield2[1]),
	    div(1.0*vfield2[0],vfield2[1]));
  test_expr(div(vfield2[0],vfield2[1]),
	    div(1.0*vfield2[0],1.0*vfield2[1]));
  test_expr(div(vfield2[0],vfield2[1]),
	    div(1.0*vfield2[0],const_cast<const array_2v&>(vfield2)[1]));
  test_expr(div(vfield2[0],vfield2[1]),
	    div(const_cast<const array_2v&>(vfield2)[0],
		const_cast<const array_2v&>(vfield2)[1]));

  // defined with BZ_ET_STENCILM. 
  test_mexpr(Jacobian3D(vfield3),
	     Jacobian3D(const_cast<const array_3v&>(vfield3)));
  test_mexpr(Jacobian3D(vfield3),
	     Jacobian3D(1.0*vfield3));
  test_mexpr(Jacobian3D(const_cast<const array_3v&>(vfield3)),
	     Jacobian3D(1.0*vfield3));

  // defined with BZ_ET_STENCILV
  test_vexpr(grad3D(field3),
	     grad3D(const_cast<const array_3&>(field3)));
  test_vexpr(grad3D(field3),
	     grad3D(1.0*field3));

  // defined with BZ_ET_STENCIL_SCA
  test_expr(div2D(vfield2),
	    div2D(const_cast<const array_2v&>(vfield2)));
  test_expr(div2D(vfield2),
	    div2D(1.0*vfield2));

  // defined with BZ_ET_STENCIL_DIFF
  test_expr(central12(field3, firstDim),
	    central12(const_cast<const array_3&>(field3), firstDim));
  test_expr(central12(field3, firstDim),
	    central12(1.0*field3, firstDim));

  result2(_bz_shrinkDomain(result2.domain(),shape(0,-1),shape(0,1))) =
    central12(where(field2>0.5,0.,1.), secondDim);
  test_expr(result2(_bz_shrinkDomain(result2.domain(),shape(0,-1),shape(0,1))), 
	    central12(where(field2>0.5,0.,1.), secondDim));
  test_expr(where(central12(field2,firstDim)>0.5,0.,1.),
	    where(central12(2*field2, firstDim)>1,0.,1.));
  test_expr(where(central12(field2, firstDim)>0.5, 
		  0.0*field2(_bz_shrinkDomain(field2.domain(),shape(-1,0),shape(1,0))),
		  0.0*field2(_bz_shrinkDomain(field2.domain(),shape(-1,0),shape(1,0)))+2.0),
	    2*where(central12(2*field2, firstDim)>1.0, 0., 1.));
  test_expr(central12(sin(1.0*field3),thirdDim), 
	    central12(1.0*sin(field3), thirdDim));
  result2 = pow(field3(Range::all(), 1, Range::all()), field2);
  test_expr(central12(result2, secondDim), 
	    central12(pow(field3(Range::all(), 1, Range::all()), 1.0*field2), secondDim));

  // defined with BZ_ET_STENCIL_MULTIDIFF
  test_expr(central12(vfield3, firstDim, secondDim),
	    central12(const_cast<const array_3v&>(vfield3), firstDim, secondDim));
  test_expr(central12(vfield3, firstDim, secondDim),
	    central12(1.0*vfield3, firstDim, secondDim));
  test_expr(where(central12(vfield2, firstDim, secondDim)>0.5, 
		  0.0*field2(_bz_shrinkDomain(field2.domain(),shape(0,-1),shape(0,1))),
		  0.0*field2(_bz_shrinkDomain(field2.domain(),shape(0,-1),shape(0,1)))+2.0),
	    2*where(central12(2*vfield2, firstDim, secondDim)>1.0, 0., 1.));

  // defined with BZ_ET_STENCIL_DIFF2
  test_expr(mixed22(field3, firstDim, secondDim),
	    mixed22(const_cast<const array_3&>(field3), firstDim, secondDim));
  test_expr(mixed22(field3, firstDim, secondDim),
	    mixed22(1.0*field3, firstDim, secondDim));
  test_expr(where(mixed22(field3, thirdDim, secondDim)>0.5, 
		  0.0*field3(_bz_shrinkDomain(field3.domain(),shape(0,-1,-1),shape(0,1,1))),
		  0.0*field3(_bz_shrinkDomain(field3.domain(),shape(0,-1,-1),shape(0,1,1)))+2.0),
	    2*where(mixed22(2*field3, thirdDim, secondDim)>1.0, 0., 1.));
    
    return 0;
}

