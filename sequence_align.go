package seq

import (
	"fmt"
)

// Alignment represents the result of aligning two sequences.
type Alignment struct {
	A []Residue // Reference
	B []Residue // Query
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
func NeedlemanWunsch(A, B []Residue, subst SubstMatrix) Alignment {
	// This implementation is taken from the "Needleman-Wunsch_algorithm"
	// Wikipedia article.
	// rows correspond to residues in A
	// cols correspond to residues in B

	// Initialization.
	var p int
	r, c := len(A)+1, len(B)+1
	matrix := make([]int, r*c)
	idx := subst.Alphabet.Index()
	sub := subst.Scores
	gapPenalty := sub[idx['-']][idx['-']]

	// Compute the matrix.
	for i := 0; i < r; i++ {
		matrix[i*c+0] = gapPenalty * i
	}
	for j := 0; j < c; j++ {
		matrix[0*c+j] = gapPenalty * j
	}

	var diag, sleft, sup int
	var subsub []int
	for i := 1; i < r; i++ {
		subsub = sub[idx[A[i-1]]]
		for j := 1; j < c; j++ {
			p = i*c + j
			diag = matrix[p-c-1] + subsub[idx[B[j-1]]]
			sup, sleft = matrix[p-c]+gapPenalty, matrix[p-1]+gapPenalty
			switch {
			case diag > sup && diag > sleft:
				matrix[p] = diag
			case sup > sleft:
				matrix[p] = sup
			default:
				matrix[p] = sleft
			}
		}
	}

	// Now trace an optimal path through the matrix starting at (r, c)
	aligned := newAlignment(max(r, c))
	i, j := r-1, c-1
	for i > 0 || j > 0 {
		p = i*c + j
		switch {
		case i > 0 && j > 0 &&
			matrix[p] == matrix[p-c-1]+sub[idx[A[i-1]]][idx[B[j-1]]]:
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
