package ioutils

import "io"

func ReadLine(r io.Reader) (string, error) {
	var ret []byte
	var tmp []byte
	tmp = make([]byte, 1)
	var err error
	err = nil

	for tmp[0] != 10 && err == nil {
		_, err = io.ReadFull(r, tmp)
		ret = append(ret, tmp[0])
	}

	return string(ret[:]), err
}

func Reverse(numbers []float64) {
	for i, j := 0, len(numbers)-1; i < j; i, j = i+1, j-1 {
		numbers[i], numbers[j] = numbers[j], numbers[i]
	}
}

func Choose(ang, ii, q []float64, test func(float64) bool) (Ang,
	I, Q []float64) {
	// body
	for i, s := range ang {
		if test(s) {
			Ang = append(Ang, s)
			I = append(I, ii[i])
			Q = append(Q, q[i])
		}
	}
	return
}
