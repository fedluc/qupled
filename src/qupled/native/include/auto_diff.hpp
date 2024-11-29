#ifndef AUTO_DIFF_HPP
#define AUTO_DIFF_HPP

// -----------------------------------------------------------------
// Classes for automatic differentiation
// -----------------------------------------------------------------

// Automatic differentiation for two dimensional functions
class AutoDiff2 {
public:

  double val;
  double dx;
  double dy;
  double dxx;
  double dyy;
  double dxy;

  // Constructors
  AutoDiff2(double val_,
	double dx_ = 0.0,
	double dy_ = 0.0,
	double dxx_ = 0.0,
	double dyy_ = 0.0,
	double dxy_ = 0.0)
    : val(val_),
      dx(dx_),
      dy(dy_),
      dxx(dxx_),
      dyy(dyy_),
      dxy(dxy_) {}
};

// addition operators
AutoDiff2 operator+(const AutoDiff2 &dual1, const AutoDiff2 &dual2);
AutoDiff2 operator+(const AutoDiff2 &dual, const double &scalar);
AutoDiff2 operator+(const double &scalar, const AutoDiff2 &dual);

// subtraction operators
AutoDiff2 operator-(const AutoDiff2 &dual1,const AutoDiff2 &dual2);
AutoDiff2 operator-(const AutoDiff2 &dual, double scalar);
AutoDiff2 operator-(const double &scalar, const AutoDiff2 &dual);

// multiplication operators
AutoDiff2 operator*(const AutoDiff2 &dual1,const AutoDiff2 &dual2);
AutoDiff2 operator*(const AutoDiff2 &dual, const double &scalar);
AutoDiff2 operator*(const double &scalar, const AutoDiff2 &dual);

// division operators
AutoDiff2 operator/(const AutoDiff2 &dual1,const AutoDiff2 &dual2);
AutoDiff2 operator/(const AutoDiff2 &dual, const double &scalar);
AutoDiff2 operator/(const double &scalar, const AutoDiff2 &dual);

// exponential function
AutoDiff2 exp(const AutoDiff2 &x);

// square root function
AutoDiff2 sqrt(const AutoDiff2 &x);

// hyperbolic tangent function
AutoDiff2 tanh(const AutoDiff2 &x);


#endif
