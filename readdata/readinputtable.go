package readdata

import (
	"fmt"
	"github.com/kshmirko/phase_cleaner/ioutils"
	"log"
	"os"
)

func ReadMeasTable(fname string) (Ang, I, Q []float64, err error) {
	// body
	fin, err := os.Open(fname)
	if err != nil {
		log.Fatal(err)
	}

	line, err := ioutils.ReadLine(fin)
	if err != nil {
		log.Fatal(err)
	}
	var tmpAng, tmpI, tmpQ, Irrad float64
	fmt.Sscanf(line, "%f %f %f %f",
		&tmpAng, &tmpQ, &tmpI, &Irrad)
	Ang = append(Ang, tmpAng)
	I = append(I, tmpI/Irrad)
	Q = append(Q, tmpQ/Irrad)

	for {
		line, err = ioutils.ReadLine(fin)
		if err != nil {
			break
		}

		fmt.Sscanf(line, "%f %f %f",
			&tmpAng, &tmpQ, &tmpI)

		Ang = append(Ang, tmpAng)
		I = append(I, tmpI/Irrad)
		Q = append(Q, tmpQ/Irrad)
	}

	return Ang, I, Q, nil

}

// Read csv sheet from excel
