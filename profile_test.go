package seq

import (
	"fmt"
	"math"
	"reflect"
	"testing"
)

var pf = fmt.Printf

var shortAlpha = NewAlphabet('A', 'B', 'C')

func o(freq, coltotal, nullfreq, nulltot int) Prob {
	if freq == 0 || nullfreq == 0 {
		return MinProb
	}
	num := float64(freq) / float64(coltotal)
	den := float64(nullfreq) / float64(nulltot)
	return -Prob(math.Log(num / den))
}

func newEProbs(alpha Alphabet, probs []map[Residue]Prob) []EProbs {
	eps := make([]EProbs, 0, len(probs))
	for _, m := range probs {
		ep := NewEProbs(alpha)
		for r, p := range m {
			ep.Set(r, p)
		}
		eps = append(eps, ep)
	}
	return eps
}

var tests = []struct {
	pexpected *Profile
	fexpected *FrequencyProfile
	null      *FrequencyProfile
	seqs      []Sequence
}{
	{
		pexpected: &Profile{
			Alphabet: shortAlpha,
			Emissions: newEProbs(shortAlpha, []map[Residue]Prob{
				{'A': o(3, 3, 3, 9), 'B': o(0, 3, 3, 9), 'C': o(0, 3, 3, 9)},
				{'A': o(0, 3, 3, 9), 'B': o(3, 3, 3, 9), 'C': o(0, 3, 3, 9)},
				{'A': o(0, 3, 3, 9), 'B': o(0, 3, 3, 9), 'C': o(3, 3, 3, 9)},
			}),
		},
		fexpected: &FrequencyProfile{
			Alphabet: shortAlpha,
			Freqs: []map[Residue]int{
				{'A': 3, 'B': 0, 'C': 0},
				{'A': 0, 'B': 3, 'C': 0},
				{'A': 0, 'B': 0, 'C': 3},
			},
		},
		null: &FrequencyProfile{
			Alphabet: shortAlpha,
			Freqs:    []map[Residue]int{{'A': 3, 'B': 3, 'C': 3}},
		},
		seqs: []Sequence{
			{"1", strr("ABC")},
			{"2", strr("ABC")},
			{"3", strr("ABC")},
		},
	},
	{
		pexpected: &Profile{
			Alphabet: shortAlpha,
			Emissions: newEProbs(shortAlpha, []map[Residue]Prob{
				{'A': o(1, 3, 3, 9), 'B': o(2, 3, 3, 9), 'C': o(0, 3, 3, 9)},
				{'A': o(0, 3, 3, 9), 'B': o(1, 3, 3, 9), 'C': o(2, 3, 3, 9)},
				{'A': o(2, 3, 3, 9), 'B': o(0, 3, 3, 9), 'C': o(1, 3, 3, 9)},
			}),
		},
		fexpected: &FrequencyProfile{
			Alphabet: shortAlpha,
			Freqs: []map[Residue]int{
				{'A': 1, 'B': 2, 'C': 0},
				{'A': 0, 'B': 1, 'C': 2},
				{'A': 2, 'B': 0, 'C': 1},
			},
		},
		null: &FrequencyProfile{
			Alphabet: shortAlpha,
			Freqs:    []map[Residue]int{{'A': 3, 'B': 3, 'C': 3}},
		},
		seqs: []Sequence{
			{"1", strr("ABC")},
			{"2", strr("BCA")},
			{"3", strr("BCA")},
		},
	},
	{
		pexpected: &Profile{
			Alphabet: shortAlpha,
			Emissions: newEProbs(shortAlpha, []map[Residue]Prob{
				{'A': o(0, 3, 1, 9), 'B': o(2, 3, 2, 9), 'C': o(1, 3, 6, 9)},
				{'A': o(0, 3, 1, 9), 'B': o(0, 3, 2, 9), 'C': o(3, 3, 6, 9)},
				{'A': o(1, 3, 1, 9), 'B': o(0, 3, 2, 9), 'C': o(2, 3, 6, 9)},
			}),
		},
		fexpected: &FrequencyProfile{
			Alphabet: shortAlpha,
			Freqs: []map[Residue]int{
				{'A': 0, 'B': 2, 'C': 1},
				{'A': 0, 'B': 0, 'C': 3},
				{'A': 1, 'B': 0, 'C': 2},
			},
		},
		null: &FrequencyProfile{
			Alphabet: shortAlpha,
			Freqs:    []map[Residue]int{{'A': 1, 'B': 2, 'C': 6}},
		},
		seqs: []Sequence{
			{"1", strr("BCC")},
			{"2", strr("BCC")},
			{"3", strr("CCA")},
		},
	},
}

func TestProfile(t *testing.T) {
	for i, test := range tests {
		columns := len(test.seqs[0].Residues)
		fprof := NewFrequencyProfileAlphabet(columns, shortAlpha)
		for _, seq := range test.seqs {
			fprof.Add(seq)
		}
		prof := fprof.Profile(test.null)
		if !reflect.DeepEqual(test.pexpected, prof) {
			t.Logf("Profile test %d failed.", i)
			t.Logf("Expected\n%v\nbut got\n%v\n", test.pexpected, prof)
			t.Fail()
		}
	}
}

func TestFrequencyProfile(t *testing.T) {
	for i, test := range tests {
		columns := len(test.seqs[0].Residues)
		prof := NewFrequencyProfileAlphabet(columns, shortAlpha)
		for _, seq := range test.seqs {
			prof.Add(seq)
		}
		if !reflect.DeepEqual(test.fexpected, prof) {
			t.Logf("Freq test %d failed.", i)
			t.Logf("Expected\n%v\nbut got\n%v\n", test.fexpected, prof)
			t.Fail()
		}
	}
}

func strr(s string) []Residue {
	bs := []byte(s)
	rs := make([]Residue, len(bs))
	for i := range bs {
		rs[i] = Residue(bs[i])
	}
	return rs
}
