package mie

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"math"
)

var RayEv *mat.Dense

func init() {
	fmt.Println("Initialize Rayleigh module")
	r := NewPowerLaw(0.0001, 0.0002, 10, 1.5-0i, -4.0, 0.75)
	_, _, _, RayEv = r.Mie("phase.mie")

}

func AotRayleigh(wl float64) float64 {
	return 0.008569 / math.Pow(wl, 4.0) *
		(1.0 + 0.0113/math.Pow(wl, 2.0) + 0.00013/math.Pow(wl, 4.0))
}
