# pytspack
 Python Wrapper around the Fortran TSPACK.

 The primary purpose of TSPACK is to construct a smooth
 function which interpolates a discrete set of data points.
 The function may be required to have either one or two con-
 tinuous derivatives, and, in the C-2 case, several options
 are provided for selecting end conditions.  If the accuracy
 of the data does not warrant interpolation, a smoothing func-
 tion (which does not pass through the data points) may be
 constructed instead.  The fitting method is designed to avoid
 extraneous inflection points (associated with rapidly varying
 data values) and preserve local shape properties of the data
 (monotonicity and convexity), or to satisfy the more general
 constraints of bounds on function values or first derivatives.
 The package also provides a parametric representation for con-
 structing general planar curves and space curves.

 The fitting function h(x) (or each component h(t) in the
 case of a parametric curve) is defined locally, on each
 interval associated with a pair of adjacent abscissae (knots),
 by its values and first derivatives at the endpoints of the
 interval, along with a nonnegative tension factor SIGMA
 associated with the interval (h is a Hermite interpolatory
 tension spline).  With SIGMA = 0, h is the cubic function
 defined by the endpoint values and derivatives, and, as SIGMA
 increases, h approaches the linear interpolant of the endpoint
 values.  Since the linear interpolant preserves positivity,
 monotonicity, and convexity of the data, h can be forced to
 preserve these properties by choosing SIGMA sufficiently
 large.  Also, since SIGMA varies with intervals, no more
 tension than necessary is used in each interval, resulting in
 a better fit and greater efficiency than is achieved with a
 single constant tension factor.
