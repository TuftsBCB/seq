package seq

import (
	"encoding/json"
)

// Alphabet corresponds to a set of residues, in a particular order, that
// capture all possible residues of a particular set of sequences. For example,
// this is used in frequency profiles and HMMs to specify which amino acids
// are represented in the probabilistic model.
//
// In most cases, the ordering of the alphabet is significant. For example,
// the indices of an alphabet may be in correspondence with the indices of
// a column in a frequency profile.
type Alphabet []Residue

// NewAlphabet creates an alphabet from the residues given.
func NewAlphabet(residues ...Residue) Alphabet {
	return Alphabet(residues)
}

func (a Alphabet) Len() int {
	return len(a)
}

// Index returns a constant-time mapping from ASCII to residue index in the
// alphabet. This depends on all residues in the alphabet being ASCII
// characters.
func (a Alphabet) Index() [256]int {
	var index [256]int
	for i, r := range a {
		index[r] = i
	}
	return index
}

// Equals returns true if and only if a1 == a2.
func (a1 Alphabet) Equals(a2 Alphabet) bool {
	if len(a1) != len(a2) {
		return false
	}
	for i, residue := range a1 {
		if residue != a2[i] {
			return false
		}
	}
	return true
}

func (a Alphabet) String() string {
	bs := make([]byte, len(a))
	for i, residue := range a {
		bs[i] = byte(residue)
	}
	return string(bs)
}

func (a *Alphabet) MarshalJSON() ([]byte, error) {
	return json.Marshal(a.String())
}

func (a *Alphabet) UnmarshalJSON(bs []byte) error {
	var str string
	if err := json.Unmarshal(bs, &str); err != nil {
		return err
	}

	*a = make(Alphabet, len(str))
	for i := 0; i < len(str); i++ {
		(*a)[i] = Residue(str[i])
	}
	return nil
}

// The default alphabet that corresponds to the BLOSUM62 matrix included
// in this package.
var AlphaBlosum62 = NewAlphabet(
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
	'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', '-',
)

// The default alphabet for DNA sequences.
var AlphaDNA = NewAlphabet(
	'A', 'C', 'G', 'T', 'N', '-',
)

// The default alphabet for RNA sequences.
var AlphaRNA = NewAlphabet(
	'A', 'C', 'G', 'U', 'N', '-',
)
