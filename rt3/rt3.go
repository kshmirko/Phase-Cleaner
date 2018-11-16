package rt3

import (
	"errors"
	"fmt"
	"github.com/kshmirko/phase_cleaner/ioutils"
	interp "github.com/kshmirko/phase_cleaner/math"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
)

type SourceType int
type DeltaScalingType string
type SurfaceType string
type RadUnitsType string
type PolarizationType string
type QuadratureType string

const (
	NoSource      SourceType = 0
	SolarSource   SourceType = 1
	ThermalSource SourceType = 2
	BothSources   SourceType = 3

	DeltaM_ON  DeltaScalingType = "Y"
	DeltaM_OFF DeltaScalingType = "N"

	STLambertian SurfaceType = "L"
	STFresnel    SurfaceType = "F"

	RUWatt           RadUnitsType = "W"
	RUBrightnessTemp RadUnitsType = "T"
	RURayleighJeans  RadUnitsType = "R"

	PTStokes        PolarizationType = "IQ"
	PTVH            PolarizationType = "VH"
	GaussQuadrature QuadratureType   = "G"
	//RT3EXEC         string           = "./rt3-code"
)

type RT3Obj struct {
	NumStokes            int
	NumQuadAngles        int
	QuadType             QuadratureType
	OrderAzimExp         int
	Layerfile            string
	DeltaM               DeltaScalingType
	Source               SourceType
	DirectFlux           float64
	DirectFluxDirection  float64
	GroundTemperature    float64
	GroundType           SurfaceType
	GroundAlbedo         float64
	SkyTemperature       float64
	Wavelengths          float64
	OutputRadUnits       RadUnitsType
	PolType              PolarizationType
	NumberOutputLayers   int
	OutputLayerIndex     int
	NumberOutputAzimuths int
	OutputFile           string
}

type RT3Res struct {
	Ang, I, Q []float64
}

func NewRt3(dir_flux, sol_zenith, galbedo, wavelen float64) *RT3Obj {
	return &RT3Obj{
		NumStokes:            2,
		NumQuadAngles:        32,
		QuadType:             GaussQuadrature,
		OrderAzimExp:         6,
		Layerfile:            "atmos.lay",
		DeltaM:               DeltaM_OFF,
		Source:               SolarSource,
		DirectFlux:           dir_flux,
		DirectFluxDirection:  sol_zenith,
		GroundTemperature:    0.0,
		GroundType:           STLambertian,
		GroundAlbedo:         galbedo,
		SkyTemperature:       0.0,
		Wavelengths:          wavelen,
		OutputRadUnits:       RUWatt,
		PolType:              PTStokes,
		NumberOutputLayers:   1,
		OutputLayerIndex:     2,
		NumberOutputAzimuths: 2,
		OutputFile:           "rt3.out",
	}
}

func (rt *RT3Obj) Run() *RT3Res {

	if ft, err := os.Stat(RT3EXEC); os.IsNotExist(err) || ft.IsDir() {
		log.Fatal(
			errors.New("There is no executable file named " + RT3EXEC),
		)
	}

	cmd := exec.Command(RT3EXEC)
	stdin, err := cmd.StdinPipe()

	if err != nil {
		log.Fatal(err)
	}

	stdout, err := cmd.StdoutPipe()
	if err != nil {
		log.Fatal(err)
	}

	err = cmd.Start()
	if err != nil {
		log.Fatal(err)
	}

	defer stdin.Close()
	defer stdout.Close()

	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.NumStokes))
	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.NumQuadAngles))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.QuadType))
	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.OrderAzimExp))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.Layerfile))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.DeltaM))
	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.Source))
	io.WriteString(stdin, fmt.Sprintf("%f\n", rt.DirectFlux))
	io.WriteString(stdin, fmt.Sprintf("%f\n", rt.DirectFluxDirection))
	io.WriteString(stdin, fmt.Sprintf("%f\n", rt.GroundAlbedo))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.GroundType))
	io.WriteString(stdin, fmt.Sprintf("%f\n", rt.GroundAlbedo))
	io.WriteString(stdin, fmt.Sprintf("%f\n", rt.SkyTemperature))
	io.WriteString(stdin, fmt.Sprintf("%f\n", rt.Wavelengths))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.OutputRadUnits))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.PolType))
	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.NumberOutputLayers))
	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.OutputLayerIndex))
	io.WriteString(stdin, fmt.Sprintf("%d\n", rt.NumberOutputAzimuths))
	io.WriteString(stdin, fmt.Sprintf("%s\n", rt.OutputFile))

	err = cmd.Wait()
	if err != nil {
		log.Fatal(err)
	}

	f, err := os.Open(rt.OutputFile)
	if err != nil {
		log.Fatal(err)
	} else {
		//log.Printf("Model evaluation successfully completed!\n")
	}

	defer f.Close()

	//Skip header from output file
	for i := 0; i < 11; i++ {
		_, _ = ioutils.ReadLine(f)
	}

	NLines := rt.NumQuadAngles * rt.NumberOutputAzimuths
	phi := make([]float64, NLines)
	mu := make([]float64, NLines)
	I := make([]float64, NLines)
	Q := make([]float64, NLines)
	ang := make([]float64, NLines)
	deg2rad := math.Pi / 180.0
	var z float64

	j := 0
	for i := 0; i < NLines*2; i++ {
		fmt.Fscanf(f, "%f %f %f %f %f\n", &z, &phi[j], &mu[j],
			&I[j], &Q[j])

		if mu[j] > 0.0 {

			ang[j] = math.Cos(phi[j]*deg2rad) *
				math.Acos(mu[j]) / deg2rad

			j++
		}
	}

	for i, j := rt.NumQuadAngles,
		2*rt.NumQuadAngles-1; i < j; i, j = i+1, j-1 {

		ang[i], ang[j] = ang[j], ang[i]
		I[i], I[j] = I[j], I[i]
		Q[i], Q[j] = Q[j], Q[i]

	}

	for i, _ := range ang {
		ang[i] = rt.DirectFluxDirection - ang[i]
	}
	// Function to filter negative values
	mytest := func(in float64) bool {
		return in > 0
	}

	ang, I, Q = ioutils.Choose(ang, I, Q, mytest)
	interp.AnglesToCos(ang)
	//ncomment if working with cosines
	//ioutils.Reverse(I)
	//ioutils.Reverse(Q)
	//ioutils.Reverse(ang)

	// Returns
	// Ang cosine of scattering angle. Ang[0] is close to 1.0 ang degrease towards -1.0, so if we whant to use returned data with other codes, we must reverse all arrays in Rt3Res
	return &RT3Res{
		I:   I,
		Q:   Q,
		Ang: ang,
	}
}

func (ret *RT3Res) Reverse() {
	// inplace reverse
	ioutils.Reverse(ret.I)
	ioutils.Reverse(ret.Q)
	ioutils.Reverse(ret.Ang)
}
