package actions

import (
	"fmt"
	"github.com/kshmirko/phase_cleaner/ioutils"
	mathutils "github.com/kshmirko/phase_cleaner/math"
	"github.com/kshmirko/phase_cleaner/mie"
	"github.com/kshmirko/phase_cleaner/readdata"
	"github.com/kshmirko/phase_cleaner/rt3"
	"github.com/urfave/cli"
	"gonum.org/v1/gonum/mat"
	"log"
	"math"
)

func InverseCalcAction(ctx *cli.Context) error {
	galbedo := ctx.Float64("galbedo")
	sol_zenith := ctx.Float64("sza")
	npts := ctx.Int("npts")
	r0 := ctx.Float64("r0")
	r1 := ctx.Float64("r1")
	gamma := ctx.Float64("gamma")
	wl := ctx.Float64("wavelen")
	midx := complex(ctx.Float64("mre"), ctx.Float64("mim"))
	taua := ctx.Float64("taua")

	Inverse(mathutils.MAX_LEG, r0, r1, npts, midx, gamma, wl, taua, sol_zenith, galbedo)
	return nil
}

func Inverse(MAXN int, R0, R1 float64, Npts int, Midx complex128,
	Gamma, Wavelen, TauA, SolarZenithDeg, GroundAlbedo float64) {

	TauM := mie.AotRayleigh(Wavelen)
	Omegas := []float64{0.5, 0.625, 0.750, 0.875,
		1.0, 1.125, 1.250}
	NOmegas := len(Omegas)

	deg2rad := math.Pi / 180.0
	mu := math.Cos(SolarZenithDeg * deg2rad)
	transmittanse := math.Exp(-(TauA + TauM) / mu)

	// Входные данные
	Ameas, Imeas, Qmeas, _ := readdata.ReadMeasTable("meas.tmp")
	mathutils.AnglesToCos(Ameas)
	ioutils.Reverse(Ameas)
	ioutils.Reverse(Imeas)
	ioutils.Reverse(Qmeas)

	//r := rand.New(rand.NewSource(99))
	for i := 0; i < len(Ameas); i++ {
		Imeas[i] = Imeas[i] * transmittanse //+ 0.005*r.Float64()
		Qmeas[i] = Qmeas[i] * transmittanse
	}

	//Prepare cosine grid

	Xi, Wi := mathutils.GaulegUnit(MAXN)

	// Prepare initial guess for size distribution function
	m := mie.NewPowerLaw(R0, R1, Npts, Midx, Gamma, Wavelen)
	_, _, _, EvA := m.Mie("phase.dat")

	//prepare phase functions
	pA_coefs := mat.Col(nil, 0, EvA)
	qA_coefs := mat.Col(nil, 1, EvA)
	pR_coefs := mat.Col(nil, 0, mie.RayEv)
	qR_coefs := mat.Col(nil, 1, mie.RayEv)

	pA := mathutils.Legval_A(Xi, pA_coefs)
	qA := mathutils.Legval_A(Xi, qA_coefs)
	pR := mathutils.Legval_A(Xi, pR_coefs)
	qR := mathutils.Legval_A(Xi, qR_coefs)

	pA1 := make([]float64, MAXN)
	qA1 := make([]float64, MAXN)
	DL := make([]float64, MAXN)

	//We bring our measurements to a single cosine grid.
	Lp := mathutils.Pwl_value_1d(Xi, Ameas, Imeas)
	Lq := mathutils.Pwl_value_1d(Xi, Ameas, Qmeas)
	_ = Lp
	_ = Lq
	_ = Xi
	_ = Wi
	log.Printf(
		"We have to iterate over single scattering albedos %d times\n",
		NOmegas)

	// Prepare RtCode Evaluators

	//rad2deg := 1.0 / deg2rad
	//mu := math.Cos(SolarZenithDeg * deg2rad)
	//transmittanse := math.Exp(-(TauA + TauM) / mu)
	log.Printf("Transmittance = %.3f\n", transmittanse)
	flux := 1.0 / mu

	rt3app0 := rt3.NewRt3(flux, SolarZenithDeg, 0.0, Wavelen)
	rt3app1 := rt3.NewRt3(flux, SolarZenithDeg, GroundAlbedo, Wavelen)

	//interp.
	// Main Loop
	var A, B float64
	norm := make([]float64, NOmegas)

	for idx, ssa_i := range Omegas {
		_, ssat, _ := rt3.PrepareScatFile(EvA, TauA, ssa_i,
			Wavelen)

		log.Printf(
			"Iteration %1d, aerosol ssa=%5.3f, total ssa=%5.3f\n",
			idx, ssa_i, ssat)

		res0 := rt3app0.Run()
		res1 := rt3app1.Run()
		res1.Reverse()
		res0.Reverse()

		L0_model := mathutils.Pwl_value_1d(Xi, res0.Ang, res0.I)
		L1_model := mathutils.Pwl_value_1d(Xi, res1.Ang, res1.I)

		for i := 0; i < MAXN; i++ {
			DL[i] = L1_model[i] - L0_model[i]
			pA1[i] = pA[i]
		}

		for i := 0; i < MAXN; i++ {
			if Lp[i] != 0.0 {
				A = ((Lp[i] - DL[i]) / L0_model[i]) * pA[i]
				B = ((Lp[i] - L1_model[i]) / L0_model[i]) *
					TauM * pR[i] / (ssa_i * TauA)
				pA1[i] = A + B
			}
		}

		int_phase0 := 0.0
		int_phase1 := 0.0
		for i := 0; i < MAXN; i++ {
			int_phase0 = int_phase0 + pA[i]*Wi[i]
			int_phase1 = int_phase1 + pA1[i]*Wi[i]

		}
		log.Printf("Norm of pA = %7.3f/%7.3f", int_phase1, int_phase0)
		norm[idx] = int_phase1

		fmt.Println()
	}

	// Find ssa where norm is 2.0
	ioutils.Reverse(Omegas)
	ioutils.Reverse(norm)

	ssa_r := mathutils.Pwl_value_1d([]float64{2.0}, norm, Omegas)[0]
	fmt.Printf("ssa_r=%f\n", ssa_r)

	if (ssa_r == 0) || (ssa_r > 2.0) {
		log.Fatal("Reconstructed omega is incorrect, please change initial" +
			" parameters")
	} else {
		_, ssat, _ := rt3.PrepareScatFile(EvA, TauA, ssa_r,
			Wavelen)

		log.Printf(
			"Final interation, aerosol ssa=%5.3f, total ssa=%5.3f\n",
			ssa_r, ssat)

		res0 := rt3app0.Run()
		res1 := rt3app1.Run()
		res1.Reverse()
		res0.Reverse()

		L0_model := mathutils.Pwl_value_1d(Xi, res0.Ang, res0.I)
		L1_model := mathutils.Pwl_value_1d(Xi, res1.Ang, res1.I)
		Q0_model := mathutils.Pwl_value_1d(Xi, res0.Ang, res0.Q)
		//Q1_model := mathutils.Pwl_value_1d(Xi, res1.Ang, res1.Q)

		for i := 0; i < MAXN; i++ {
			DL[i] = L1_model[i] - L0_model[i]
			pA1[i] = pA[i]
			qA1[i] = qA[i]
		}

		for i := 0; i < MAXN; i++ {
			if Lp[i] != 0.0 {
				A = ((Lp[i] - DL[i]) / L0_model[i]) * pA[i]
				B = ((Lp[i] - L1_model[i]) / L0_model[i]) *
					TauM * pR[i] / (ssa_r * TauA)
				pA1[i] = A + B

				A = (Lq[i] / Q0_model[i]) * qA[i]
				B = ((Lq[i] - Q0_model[i]) / Q0_model[i]) *
					TauM * qR[i] / (ssa_r * TauA)
				qA1[i] = A + B
			}
		}

		int_phase0 := 0.0
		int_phase1 := 0.0
		for i := 0; i < MAXN; i++ {
			int_phase0 = int_phase0 + pA[i]*Wi[i]
			int_phase1 = int_phase1 + pA1[i]*Wi[i]

		}
		log.Printf("Norm of pA = %7.3f/%7.3f", int_phase1, int_phase0)
	}

	for i := 0; i < MAXN; i++ {
		pA1[i] = (TauA*ssa_r*pA1[i] + TauM*pR[i]) / (TauM + TauA*ssa_r)
		qA1[i] = (TauA*ssa_r*qA1[i] + TauM*qR[i]) / (TauM + TauA*ssa_r)
		fmt.Printf("%10.3f%12.6f%12.6f%12.6f%12.6f\n",
			Xi[i], pA1[i], qA1[i], pA[i], qA[i])
	}

}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%12.4e\n", fa)
}
