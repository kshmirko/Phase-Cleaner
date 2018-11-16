package mathutils

import (
	"math"
)

const (
	MAX_LEG int = 64
)

func AnglesToCos(ang []float64) {
	deg2rad := math.Pi / 180.0
	for i, mu_i := range ang {
		ang[i] = math.Cos(mu_i * deg2rad)
	}
	return
}

func Gauleg(x1, x2 float64, npts int) (xi, wi []float64) {
	// body
	const EPS float64 = 3e-11
	var m, j, i int
	var z1, z, xm, xl, pp, p3, p2, p1 float64
	xi = make([]float64, npts)
	wi = make([]float64, npts)
	m = (npts + 1) / 2
	xm = 0.5 * (x2 + x1)
	xl = 0.5 * (x2 - x1)
	for i = 1; i <= m; i++ {
		z = math.Cos(math.Pi * (float64(i) - 0.25) /
			(float64(npts) + 0.5))
		for {
			p1 = 1.0
			p2 = 0.0
			for j = 1; j <= npts; j++ {
				p3 = p2
				p2 = p1
				p1 = ((2.0*float64(j)-1.0)*z*p2 - (float64(j)-1.0)*p3) /
					float64(j)
			}
			pp = float64(npts) * (z*p1 - p2) / (z*z - 1.0)
			z1 = z
			z = z1 - p1/pp
			if math.Abs(z-z1) < EPS {
				break
			}
		}
		xi[i-1] = xm - xl*z
		xi[(npts+1-i)-1] = xm + xl*z
		wi[i-1] = 2.0 * xl / ((1.0 - z*z) * pp * pp)
		wi[(npts+1-i)-1] = wi[i-1]
	}
	return
}

func GaulegUnit(npts int) (xi, wi []float64) {
	return Gauleg(-1, 1, npts)
}

func Pwl_value_1d(xi, xd, yd []float64) (yi []float64) {

	/******************************************************************************/
	/*
	   Purpose:

	     PWL_VALUE_1D evaluates the piecewise linear interpolant.

	   Discussion:

	     The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
	     linear function which interpolates the data (XD(I),YD(I)) for I = 1
	     to ND.

	   Licensing:

	     This code is distributed under the GNU LGPL license.

	   Modified:

	     22 September 2012

	   Author:

	     John Burkardt

	   Parameters:

	     Input, int ND, the number of data points.
	     ND must be at least 1.

	     Input, double XD[ND], the data points.

	     Input, double YD[ND], the data values.

	     Input, int NI, the number of interpolation points.

	     Input, double XI[NI], the interpolation points.

	     Output, double PWL_VALUE_1D[NI], the interpolated values.
	*/
	var i, k int
	var t float64
	ni := len(xi)
	nd := len(xd)
	yi = make([]float64, ni)

	if nd == 1 {
		for i = 0; i < ni; i++ {
			yi[i] = yd[0]
		}
		return
	}

	for i = 0; i < ni; i++ {
		if xi[i] <= xd[0] {
			yi[i] = 0.0
		} else if xd[nd-1] <= xi[i] {
			yi[i] = 0.0
		} else {
			for k = 1; k < nd; k++ {
				if xd[k-1] <= xi[i] && xi[i] <= xd[k] {
					t = (xi[i] - xd[k-1]) / (xd[k] - xd[k-1])
					yi[i] = (1.0-t)*yd[k-1] + t*yd[k]
					break
				}
			}
		}
	}
	return
}

func Legval(x float64, coeff []float64) (legval float64) {
	n := len(coeff) - 1

	P0 := 1.0
	P1 := x
	P2 := 0.0

	legval = coeff[0]*P0 + coeff[1]*P1

	for i := 2; i <= n; i++ {
		P2 = float64(2*i-1)/float64(i)*x*P1 - float64(i-1)/float64(i)*P0
		legval = legval + coeff[i]*P2
		P0 = P1
		P1 = P2
	}

	return
}

func Legval_A(x, coeff []float64) []float64 {
	// body
	m := len(x)
	y := make([]float64, m)

	for i := 0; i < m; i++ {
		y[i] = Legval(x[i], coeff)
	}
	return y
}
