package mie

import (
	"errors"
	"fmt"
	"github.com/kshmirko/phase_cleaner/ioutils"
	"gonum.org/v1/gonum/mat"
	"io"
	"log"
	"os"
	"os/exec"
)

type PowerLawDistr struct {
	R0, R1  float64
	Nbins   int
	Midx    complex128
	Gamma   float64
	Wavelen float64
}

func NewPowerLaw(r0, r1 float64, nbins int, midx complex128,
	gamma, wl float64) *PowerLawDistr {
	return &PowerLawDistr{
		R0:      r0,
		R1:      r1,
		Nbins:   nbins,
		Midx:    midx,
		Gamma:   gamma,
		Wavelen: wl,
	}
}

func (distr *PowerLawDistr) Mie(fname string) (float64, float64, float64, *mat.Dense) {
	if ft, err := os.Stat(MIEEXE); os.IsNotExist(err) || ft.IsDir() {
		log.Fatal(
			errors.New("There is no executable file named " + MIEEXE),
		)
	}

	cmd := exec.Command(MIEEXE)
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

	io.WriteString(stdin,
		fmt.Sprintf("%f\n%f\n%d\n%f\n(%f,%f)\n%f\n",
			distr.R0,
			distr.R1,
			distr.Nbins,
			distr.Gamma,
			real(distr.Midx), imag(distr.Midx),
			distr.Wavelen,
		),
	)

	// Read stdout
	ext, sca, ssa, mat := distr.read_stdout(fname, stdout)

	//wait for command finished
	err = cmd.Wait()

	if err != nil {
		log.Fatal(err)
	}

	return ext, sca, ssa, mat
}

func (distr *PowerLawDistr) read_stdout(fname string, stdout io.Reader) (float64,
	float64, float64, *mat.Dense) {

	fout, err := os.Create(fname)
	if err != nil {
		log.Fatal(err)
	}

	defer fout.Close()

	// Read the output
	for i := 0; i < 6; i++ {
		_, _ = ioutils.ReadLine(stdout)
	}

	var ext, sca, ssa float64
	var nleg int

	tmp, err := ioutils.ReadLine(stdout)
	fout.WriteString(tmp)
	fmt.Sscanf(tmp, "%f", &ext)
	tmp, err = ioutils.ReadLine(stdout)
	fout.WriteString(tmp)
	fmt.Sscanf(tmp, "%f", &sca)
	tmp, err = ioutils.ReadLine(stdout)
	fout.WriteString(tmp)
	fmt.Sscanf(tmp, "%f", &ssa)

	tmp, err = ioutils.ReadLine(stdout)
	fout.WriteString(tmp)
	fmt.Sscanf(tmp, "%d", &nleg)

	mat := mat.NewDense(nleg+1, 6, nil)
	row := make([]float64, 6)

	// read evans matrix
	for i := 0; i <= nleg; i++ {

		tmp, err = ioutils.ReadLine(stdout)
		if err != nil {
			log.Fatal(err)
		}
		fout.WriteString(tmp)

		var idx int

		fmt.Sscanf(tmp, "%d %f %f %f %f %f %f", &idx, &row[0],
			&row[1], &row[2], &row[3], &row[4], &row[5])

		mat.SetRow(idx, row)
	}

	return ext, sca, ssa, mat
}
