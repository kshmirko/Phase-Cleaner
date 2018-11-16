package miecode

/*
#cgo CFLAGS: -Os
#cgo LDFLAGS: -L. -lmie
#include "libmie.h"
*/
import (
	"C"
)
import (
	_ "fmt"
	mathutils "github.com/kshmirko/phase_cleaner/math"
	"gonum.org/v1/gonum/mat"
	"math"
	"math/cmplx"
	"unsafe"
)

func toCArray(x []float64) []C.double {
	v := make([]C.double, len(x))
	for i := 0; i < len(x); i++ {
		v[i] = (C.double)(x[i])
	}
	return v
}

func toGoArray(x []C.double) []float64 {
	v := make([]float64, len(x))
	for i := 0; i < len(x); i++ {
		v[i] = float64(x[i])
	}
	return v
}

type MieCode int

func (mie *MieCode) EazyMie(x float64, n float64) (qsca float64, g float64) {
	C.ez_Mie(C.double(x), C.double(n),
		(*C.double)(unsafe.Pointer(&qsca)),
		(*C.double)(unsafe.Pointer(&g)))
	return qsca, g
}

func (mie *MieCode) EazyMieFull(x float64, m complex128, mu []float64) (s1 []complex128, s2 []complex128, qext, qsca, qback, g float64) {

	mu_ := toCArray(mu)
	s1real := make([]C.double, len(mu))
	s2real := make([]C.double, len(mu))
	s1imag := make([]C.double, len(mu))
	s2imag := make([]C.double, len(mu))

	arrSize := len(mu)
	C.ez_Mie_Full(C.double(x), C.double(real(m)), C.double(imag(m)),
		C.long(arrSize), (*C.double)(&mu_[0]),
		(*C.double)(unsafe.Pointer(&s1real[0])),
		(*C.double)(unsafe.Pointer(&s1imag[0])),
		(*C.double)(unsafe.Pointer(&s2real[0])),
		(*C.double)(unsafe.Pointer(&s2imag[0])),
		(*C.double)(unsafe.Pointer(&qext)),
		(*C.double)(unsafe.Pointer(&qsca)),
		(*C.double)(unsafe.Pointer(&qback)),
		(*C.double)(unsafe.Pointer(&g)))

	s1 = make([]complex128, arrSize)
	s2 = make([]complex128, arrSize)

	for i := 0; i < arrSize; i++ {
		s1[i] = complex(s1real[i], s1imag[i])
		s2[i] = complex(s2real[i], s2imag[i])
	}
	return s1, s2, qext, qsca, qback, g
}

func (mie *MieCode) MuellerMatrixMie(x float64, m complex128, mu []float64) (ret *mat.Dense, qext, qsca, qback, g float64) {

	s1, s2, qext, qsca, qback, g := mie.EazyMieFull(x, m, mu)
	numang := len(mu)
	ret = mat.NewDense(numang, 6, nil)

	row := []float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	Sqrt4Pi := complex(math.Sqrt(4*math.Pi), 0)
	for i, _ := range s1 {
		s1[i] = s1[i] * Sqrt4Pi
		s2[i] = s2[i] * Sqrt4Pi
		row[0] = 0.5 * (math.Pow(cmplx.Abs(s2[i]), 2) +
			math.Pow(cmplx.Abs(s1[i]), 2))
		row[1] = 0.5 * (math.Pow(cmplx.Abs(s2[i]), 2) -
			math.Pow(cmplx.Abs(s1[i]), 2))
		row[2] = real(cmplx.Conj(s2[i]) * s1[i])
		row[3] = imag(s2[i] * cmplx.Conj(s1[i]))
		row[4] = row[0]
		row[5] = row[2]

		ret.SetRow(i, row)
	}

	return
}

func (mie *MieCode) P_polynomial_value(n int, x []float64) *mat.Dense {

	m := len(x)

	ret := mat.NewDense(n, m, nil)

	for i := 0; i < m; i++ {
		ret.Set(0, i, 1.0)
	}

	if n < 1 {
		return ret
	}

	for i := 0; i < m; i++ {
		ret.Set(1, i, x[i])
	}

	for j := 2; j < n; j++ {
		for i := 0; i < m; i++ {
			tmp := (float64(2*j-1)*x[i]*ret.At(j-1, i) -
				float64(j-1)*ret.At(j-2, i)) / float64(j)
			ret.Set(j, i, tmp)
		}
	}

	return ret
}

func (mie *MieCode) FindLegendreCoefficients(maxdeg int, x []float64, evans *mat.Dense) *mat.Dense {

	A := mie.P_polynomial_value(maxdeg, x)
	ra, _ := A.Dims()
	_, ce := evans.Dims()
	ret := mat.NewDense(ra, ce, nil)
	if err := ret.Solve(A.T(), evans); err != nil {
		return nil
	}
	return ret
}

func (mie *MieCode) CalcScatteringProps(r0, r1, gamma, wavelen float64, nleg int, midx complex128, mu []float64) (evans_tot *mat.Dense, ext, sca, back, ssa float64) {

	f := func(r, gamma float64) float64 { return math.Pow(r, gamma) }
	norm := f(r1, gamma+1)/(gamma+1) - f(r0, gamma+1)/(gamma+1)
	k := 2 * math.Pi / wavelen
	//dr := (r1 - r0) / (float64(npts) - 1)

	ri, wi := mathutils.Gauleg(r0, r1, mathutils.MAX_LEG)

	var r_i, f_i, s_i, x_i float64
	var csca_i float64
	evans_tot = mat.NewDense(nleg, 6, nil)

	for i := 0; i < mathutils.MAX_LEG; i++ {
		r_i = ri[i]
		f_i = f(r_i, gamma) / norm
		s_i = math.Pi * r_i * r_i
		x_i = k * r_i
		evans, qext, qsca, qback, _ :=
			mie.MuellerMatrixMie(x_i, midx, mu)
		csca_i = qsca * s_i * f_i
		ext = ext + qext*s_i*f_i*wi[i]
		sca = sca + csca_i*wi[i]
		back = back + qback*s_i*f_i*wi[i]

		Evans := mie.FindLegendreCoefficients(nleg, mu, evans)

		Evans.Scale(csca_i*wi[i], Evans)
		evans_tot.Add(evans_tot, Evans)
	}

	evans_tot.Scale(1.0/evans_tot.At(0, 0), evans_tot)
	ext *= 1e-12
	sca *= 1e-12
	back *= 1e-12

	return
}
