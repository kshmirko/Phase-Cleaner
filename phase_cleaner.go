package main

import (
	"github.com/kshmirko/phase_cleaner/actions"
	"github.com/urfave/cli"
	"log"
	"os"
	"time"
	//"math/rand"
)

func main() {

	app := cli.NewApp()

	app.Commands = []cli.Command{
		cli.Command{
			Name:     "direct",
			Category: "Direct problem",
			Usage:    "Modeling solar radiance at the ground level",
			Flags: []cli.Flag{
				cli.Float64Flag{
					Name:  "r0, 0",
					Usage: "Smallest particle's size",
					Value: 0.1,
				},
				cli.Float64Flag{
					Name:  "r1, 1",
					Usage: "Largest particle's size",
					Value: 1.0,
				},
				cli.Float64Flag{
					Name:  "gamma, g",
					Usage: "Size distribution decay",
					Value: -4.0,
				},
				cli.Float64Flag{
					Name:  "wavelen, w",
					Usage: "Incident radiation wavelength",
					Value: 0.75,
				},
				cli.IntFlag{
					Name:  "npts, n",
					Usage: "Number of SD-bins",
					Value: 101,
				},
				cli.Float64Flag{
					Name:  "taua, t",
					Usage: "Aerosol optical depth",
					Value: 0.1,
				},
				cli.Float64Flag{
					Name:  "sza",
					Usage: "Solar zenith angle",
					Value: 65.0,
				},
				cli.Float64Flag{
					Name:  "galbedo",
					Usage: "Surface albedo",
					Value: 0.1,
				},
				cli.Float64Flag{
					Name:  "mre",
					Usage: "Real part of particle refractive index",
					Value: 1.5,
				},
				cli.Float64Flag{
					Name:  "mim",
					Usage: "Imaginary part of particle refractive index",
					Value: 0,
				},
			},
			Action: actions.DirectCalcAction,
		},
		cli.Command{
			Name:     "inverse",
			Category: "Inverse problem",
			Usage:    "Perform phase function cleaning from effects of multiple scattering and ground reflection",
			Flags: []cli.Flag{
				cli.Float64Flag{
					Name:  "r0, 0",
					Usage: "Smallest particle's size",
					Value: 0.1,
				},
				cli.Float64Flag{
					Name:  "r1, 1",
					Usage: "Largest particle's size",
					Value: 1.0,
				},
				cli.Float64Flag{
					Name:  "gamma, g",
					Usage: "Size distribution decay",
					Value: -4.0,
				},
				cli.Float64Flag{
					Name:  "wavelen, w",
					Usage: "Incident radiation wavelength",
					Value: 0.75,
				},
				cli.IntFlag{
					Name:  "npts, n",
					Usage: "Number of SD-bins",
					Value: 101,
				},
				cli.Float64Flag{
					Name:  "taua, t",
					Usage: "Aerosol optical depth",
					Value: 0.1,
				},
				cli.Float64Flag{
					Name:  "sza",
					Usage: "Solar zenith angle",
					Value: 65.0,
				},
				cli.Float64Flag{
					Name:  "galbedo",
					Usage: "Surface albedo",
					Value: 0.1,
				},
				cli.Float64Flag{
					Name:  "mre",
					Usage: "Real part of particle refractive index",
					Value: 1.5,
				},
				cli.Float64Flag{
					Name:  "mim",
					Usage: "Imaginary part of particle refractive index",
					Value: 0,
				},
			},
			Action: actions.InverseCalcAction,
		},
	}

	app.Action = nil
	app.Usage = "Perform cleaning solar radiance measured at ground level for multiple scattering and ground reflection"
	app.Version = "0.1.0"
	app.Compiled = time.Now()
	app.Authors = []cli.Author{
		cli.Author{
			Name:  "Konstantin Shmirko",
			Email: "kshmirko(at)gmail.com",
		},
	}
	app.Copyright = "(c) 2018 IACP FEB RAS"
	err := app.Run(os.Args)
	if err != nil {
		log.Fatal(err)
	}

	// Define command line keys
	// flag.IntVar(&Npts, "npts", 101, "Number of SD-bins")
	//     flag.IntVar(&Nleg, "nleg", 20, "Number Legendre coefficients")
	//     flag.Float64Var(&R0, "r0", 0.1, "Smallest particle's size")
	//     flag.Float64Var(&R1, "r1", 1.0, "Largest particle's size")
	//     flag.Float64Var(&Gamma, "gamma", -3.5, "Size distribution decay")
	//     flag.Float64Var(&Wavelen, "wavlen", 0.750, "Wavelengths of incident"+
	//         " radiation")
	//     flag.Float64Var(&TauA, "taua", 0.1, "Aerosol optical depth")
	//     flag.Float64Var(&SolarZenithDeg, "sza", 65, "Solar zenith angle")
	//     flag.Float64Var(&GroundAlbedo, "galbedo", 0.1, "Ground albedo")
	//     flag.Float64Var(&mre, "mre", 1.5, "Real part of particle refractive index")
	//     flag.Float64Var(&mim, "mim", 0.0, "Imaginary part of particle refractive"+
	//         " index")
	//     flag.Parse()

	// Test initial parameters
	// MAXN := 64
	//R0 := 0.1
	//R1 := 1.0
	//Npts := *nptsv
	// Midx := complex(mre, mim)
	//Gamma := -3.5
	//Wavelen := 0.75
	//TauA := 0.10
	//SolarZenithDeg := 65.0
	//GroundAlbedo := 0.1

	// var mie miecode.MieCode

	// xi, _ := mathutils.GaulegUnit(MAXN)
	// Etot, ext, sca, _, _ := mie.CalcScatteringProps(R0, R1, Gamma, Wavelen, Nleg, Midx, xi)

	// fmt.Printf("Ext = %12.3e\t Sca = %12.3e\nSsa = %7.3f\n", ext, sca, sca/ext)
	// matPrint(Etot)
	//do_direct_calc(GroundAlbedo, SolarZenithDeg, Npts, R0, R1, Gamma, Wavelen, Midx, TauA)

	// Inverse(MAXN, R0, R1, Npts, Midx, Gamma, Wavelen, TauA, SolarZenithDeg, GroundAlbedo)

}
