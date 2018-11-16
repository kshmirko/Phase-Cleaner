package main

import "testing"
import "github.com/kshmirko/phase_cleaner/actions"

func BenchmarkDirect(*testing.B) {
	actions.Do_direct_calc(0.1, 65.0, 101,
		0.1, 1.0, -3.5, 0.75, 1.5+0i, 0.1)
}

func BenchmarkInverse(*testing.B) {
	actions.Inverse(64, 0.1, 1.0, 101, 1.5+0i,
		-3.5, 0.75, 0.1, 65.0, 0.1)
}
