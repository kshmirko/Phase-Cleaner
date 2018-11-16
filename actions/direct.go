package actions

import (
	"fmt"
	_ "github.com/kshmirko/phase_cleaner/math"
	"github.com/kshmirko/phase_cleaner/mie"
	"github.com/kshmirko/phase_cleaner/rt3"
	"github.com/urfave/cli"
	"math"
)

func DirectCalcAction(ctx *cli.Context) error {
	galbedo := ctx.Float64("galbedo")
	sol_zenith := ctx.Float64("sza")
	npts := ctx.Int("npts")
	r0 := ctx.Float64("r0")
	r1 := ctx.Float64("r1")
	gamma := ctx.Float64("gamma")
	wl := ctx.Float64("wavelen")
	midx := complex(ctx.Float64("mre"), ctx.Float64("mim"))
	taua := ctx.Float64("taua")

	Do_direct_calc(galbedo, sol_zenith, npts, r0, r1, gamma, wl, midx, taua)
	return nil
}

func Do_direct_calc(galbedo, sol_zenith float64, npts int,
	r0, r1, gamma, wl float64, midx complex128, taua float64) {

	m := mie.NewPowerLaw(r0, r1, npts,
		midx, gamma, wl)

	taum := mie.AotRayleigh(wl)
	deg2rad := math.Pi / 180.0
	mu := math.Cos(sol_zenith * deg2rad)
	transmittanse := math.Exp(-(taua + taum) / mu)
	_, _, ssa, EvA := m.Mie("phase.dat")

	_, _, _ = rt3.PrepareScatFile(EvA, taua, ssa,
		wl)

	flux := 1.0 / mu
	rt3app1 := rt3.NewRt3(flux, sol_zenith, galbedo, wl)
	res1 := rt3app1.Run()
	//res1.Reverse()

	for i := 0; i < len(res1.I); i++ {
		res1.I[i] /= transmittanse
		res1.Q[i] /= transmittanse
	}

	//MuSol := math.Cos(deg2rad * sol_zenith)
	for i := 0; i < len(res1.I); i++ {
		fmt.Printf("%15.8f%15.8f%15.8f\n", math.Acos(res1.Ang[i])/deg2rad,
			res1.Q[i], res1.I[i])
	}
}
