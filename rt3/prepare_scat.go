package rt3

import (
	"github.com/kshmirko/phase_cleaner/mie"
	"gonum.org/v1/gonum/mat"

	"fmt"
	"log"
	"os"
)

func PrepareScatFile(EvA *mat.Dense,
	taua, ssa, wavelen float64) (float64, float64, *mat.Dense) {
	r, c := EvA.Dims()
	ret := mat.NewDense(r, c, nil)
	tmp := mat.NewDense(r, c, nil)
	taua_s := ssa * taua
	taum := mie.AotRayleigh(wavelen)

	ret.Scale(taua_s, EvA)
	tmp.Scale(taum, mie.RayEv)
	ret.Add(ret, tmp)
	ret.Scale(1.0/(taua_s+taum), ret)
	taut := taum + taua
	ssat := (taua_s + taum) / (taua + taum)

	fout, err := os.Create("scat_file")
	if err != nil {
		log.Fatal(err)
	}

	defer fout.Close()

	fout.WriteString(fmt.Sprintf("%12.4f\n", taut))
	fout.WriteString(fmt.Sprintf("%12.4f\n", taut*ssat))
	fout.WriteString(fmt.Sprintf("%12.4f\n", taut))
	fout.WriteString(fmt.Sprintf("%12d\n", r-1))

	for i := 0; i < r; i++ {
		fout.WriteString(
			fmt.Sprintf("%7d%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n",
				i, ret.At(i, 0), ret.At(i, 1), ret.At(i, 2),
				ret.At(i, 3), ret.At(i, 4), ret.At(i, 5),
			),
		)
	}

	return taut, ssat, ret
}
