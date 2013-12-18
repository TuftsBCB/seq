package seq

import (
	"fmt"
)

type Alignment struct {
	A, B []Residue
}

func newAlignment(length int) Alignment {
	return Alignment{
		A: make([]Residue, 0, length),
		B: make([]Residue, 0, length),
	}
}

// Performs the Needleman-Wunsch sequence alignment algorithm on a pair
// of sequences. A function `subst` should return alignment scores for
// pairs of residues. This package provides some functions suitable for
// this purpose, e.g., MatBlosum62, MatDNA, MatRNA, etc.
func NeedlemanWunsch(A, B []Residue, subst MatLookup) Alignment {
	// This implementation is taken from the "Needleman-Wunsch_algorithm"
	// Wikipedia article.
	// rows correspond to residues in A
	// cols correspond to residues in B

	// Initialization.
	r, c := len(A)+1, len(B)+1
	gapPenalty := subst('-', '-')
	matrix := make([]int, r*c)

	// Compute the matrix.
	for i := 0; i < r; i++ {
		matrix[i*c+0] = gapPenalty * i
	}
	for j := 0; j < c; j++ {
		matrix[0*c+j] = gapPenalty * j
	}
	for i := 1; i < r; i++ {
		for j := 1; j < c; j++ {
			p := i*c + j
			matrix[p] = max3(
				matrix[p-c-1]+subst(A[i-1], B[j-1]),
				matrix[p-c]+gapPenalty,
				matrix[p-1]+gapPenalty)
		}
	}

	// Now trace an optimal path through the matrix starting at (r, c)
	aligned := newAlignment(max(r, c))
	i, j := r-1, c-1
	for i > 0 || j > 0 {
		p := i*c + j
		switch {
		case i > 0 && j > 0 && matrix[p] == matrix[p-c-1]+subst(A[i-1], B[j-1]):
			aligned.A = append(aligned.A, A[i-1])
			aligned.B = append(aligned.B, B[j-1])
			i--
			j--
		case i > 0 && matrix[p] == matrix[p-c]+gapPenalty:
			aligned.A = append(aligned.A, A[i-1])
			aligned.B = append(aligned.B, '-')
			i--
		case j > 0 && matrix[p] == matrix[p-1]+gapPenalty:
			aligned.A = append(aligned.A, '-')
			aligned.B = append(aligned.B, B[j-1])
			j--
		default:
			panic(fmt.Sprintf("BUG in NeedlemanWunsch: No path at (%d, %d)",
				i, j))
		}
	}

	// Since we built the alignment in backwards, we must reverse the alignment.
	for i, j := 0, len(aligned.A)-1; i < j; i, j = i+1, j-1 {
		aligned.A[i], aligned.A[j] = aligned.A[j], aligned.A[i]
		aligned.B[i], aligned.B[j] = aligned.B[j], aligned.B[i]
	}
	return aligned
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func max3(a, b, c int) int {
	switch {
	case a > b && a > c:
		return a
	case b > c:
		return b
	}
	return c
}
