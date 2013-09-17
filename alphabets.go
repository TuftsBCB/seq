package seq

import (
	"encoding/json"
)

type Alphabet []Residue

func NewAlphabet(residues ...Residue) Alphabet {
	return Alphabet(residues)
}

func (a Alphabet) Len() int {
	return len(a)
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

const alpha62letters = "ABCDEFGHIKLMNPQRSTVWXYZ-"

// The default alphabet that corresponds to the BLOSUM62 matrix included
// in this package.
var AlphaBlosum62 = NewAlphabet(
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
	'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', '-',
)

const alphaDNAletters = "ACTGN-"
const alphaRNAletters = "ACUGN-"

var AlphaDNA = NewAlphabet(
	'A', 'C', 'G', 'T', 'N', '-',
)

var AlphaRNA = NewAlphabet(
	'A', 'C', 'G', 'U', 'N', '-',
)