package seq

import (
	"encoding/json"
	"fmt"
	"math"
	"strconv"
)

// HMM states in the Plan7 architecture.
const (
	Match HMMState = iota
	Deletion
	Insertion
	Begin
	End
)

type HMMState int

type HMM struct {
	// An ordered list of HMM nodes.
	Nodes []HMMNode

	// The alphabet as defined by an ordering of residues.
	// Indices in this slice correspond to indices in match/insertion emissions.
	Alphabet Alphabet

	// NULL model. (Amino acid background frequencies.)
	// HMMER hmm files don't have this, but HHsuite hhm files do.
	// In the case of HHsuite, the NULL model is used for insertion emissions
	// in every node.
	Null EProbs
}

// DynamicTable represents a dynamic programming table used in sequence
// alignment algorithms like Viterbi.
type DynamicTable struct {
	scores []Prob
	nodes  int
}

// AllocTable returns a freshly allocated dynamic programming table for use
// in HMM alogirthms like Viterbi. It is indexed by HMM state, node index and
// sequence length, in that order. The total size of the table is equal to
// (#states * (#nodes + 1) * (seqLen + 1)).
//
// The only states used are Match, Deletion and Insertion.
//
// Each value is initialized to a minimum probability.
func AllocTable(numNodes int, seqLen int) *DynamicTable {
	nodes := numNodes + 1
	t := &DynamicTable{
		scores: make([]Prob, 3*nodes*(seqLen+1)),
		nodes:  nodes,
	}
	t.reset()
	return t
}

func (t *DynamicTable) index(state HMMState, node int, obs int) int {
	return int(state) + 3*(node+t.nodes*obs)
}

func (t *DynamicTable) set(state HMMState, node int, obs int, p Prob) {
	i := t.index(state, node, obs)
	if t.scores[i].Less(p) {
		t.scores[i] = p
	}
}

func (t *DynamicTable) reset() {
	for i := 0; i < len(t.scores); i++ {
		t.scores[i] = MinProb
	}
}

// ViterbiScore returns the probability of the likeliest path through the HMM
// for the given sequence.
//
// If you're running Viterbi in a performance critical section, ViterbiScoreMem
// may be appropriate.
//
// Note that the state path is not computed.
func (hmm *HMM) ViterbiScore(seq Sequence) Prob {
	table := AllocTable(len(hmm.Nodes), seq.Len())
	return hmm.ViterbiScoreMem(seq, table)
}

// ViterbiScoreMem is the same as ViterbiScore, except it does not allocate,
// which makes it faster in performance critical sections of code. This is done
// by passing a pre-allocated dynamic programming table created by AllocTable
// function.
//
// Note that the caller must ensure that only one goroutine is calling
// ViterbiScoreMem with the same dynamic programming table.
func (hmm *HMM) ViterbiScoreMem(seq Sequence, table *DynamicTable) Prob {
	table.reset()
	table.scores[table.index(Match, 0, 0)] = Prob(0.0) // The begin node.

	var trans TProbs
	var residue Residue
	var memit, iemit, here Prob
	for node := 0; node < len(hmm.Nodes); node++ {
		for obs := 0; obs < seq.Len(); obs++ {
			trans = hmm.Nodes[node].Transitions
			residue = seq.Residues[obs]
			iemit = hmm.Nodes[node].InsEmit.Lookup(residue)
			if node+1 < len(hmm.Nodes) {
				memit = hmm.Nodes[node+1].MatEmit.Lookup(residue)
			} else {
				memit = 0.0 // Force into match state for end node.
			}

			here = table.scores[table.index(Match, node, obs)]
			table.set(Insertion, node, obs+1, here+trans.MI+iemit)
			table.set(Match, node+1, obs+1, here+trans.MM+memit)
			table.set(Deletion, node+1, obs, here+trans.MD)

			here = table.scores[table.index(Insertion, node, obs)]
			table.set(Insertion, node, obs+1, here+trans.II+iemit)
			table.set(Match, node+1, obs+1, here+trans.IM+memit)

			here = table.scores[table.index(Deletion, node, obs)]
			table.set(Match, node+1, obs+1, here+trans.DM+memit)
			table.set(Deletion, node+1, obs, here+trans.DD)
		}
	}
	return table.scores[table.index(Match, len(hmm.Nodes), seq.Len())]
}

// HMMNode represents a single node in an HMM, including the reference residue,
// the node index, insertion emissions, match emissions, transition
// probabilities.
//
// The NeffM, NeffI and NeffD aren't used, but are included since they exist
// in common HMM file formats.
type HMMNode struct {
	Residue             Residue
	NodeNum             int
	InsEmit             EProbs
	MatEmit             EProbs
	Transitions         TProbs
	NeffM, NeffI, NeffD Prob
}

// EProbs represents emission probabilities, as log-odds scores.
// The representation of EProbs should not be used by clients; it is exported
// only for convenience with marshalling APIs.
type EProbs struct {
	Offset Residue
	Probs  []Prob
}

// NewEProbs creates a new EProbs map from the given alphabet. Keys of the map
// are residues defined in the alphabet, and values are defaulted to the
// minimum probability.
func NewEProbs(alphabet Alphabet) EProbs {
	offset, max := Residue(255), Residue(0)
	for _, r := range alphabet {
		if r < offset {
			offset = r
		}
		if r > max {
			max = r
		}
	}
	probs := make([]Prob, 1+max-offset)
	for i := range probs {
		probs[i] = MinProb
	}
	return EProbs{offset, probs}
}

// Returns the emission probability for a particular residue.
func (ep EProbs) Lookup(r Residue) Prob {
	i := r - ep.Offset
	if i < 0 || int(i) >= len(ep.Probs) {
		return MinProb
	}
	return ep.Probs[r-ep.Offset]
}

// Set sets the emission probability of the given residue.
func (ep *EProbs) Set(r Residue, p Prob) {
	ep.Probs[r-ep.Offset] = p
}

// TProbs represents transition probabilities, as log-odds scores.
// Note that ID and DI are omitted (Plan7).
type TProbs struct {
	MM, MI, MD, IM, II, DM, DD Prob
}

// Prob represents a transition or emission probability.
type Prob float64

var invalidProb = Prob(math.NaN())

// The value representing a minimum emission/transition probability.
// Remember, max in negative log space is minimum probability.
var MinProb = Prob(math.MaxFloat64)

// NewProb creates a new probability value from a string (usually read from
// an hmm or hhm file). If the string is equivalent to the special value "*",
// then the probability returned is guaranteed to be minimal. Otherwise, the
// string is parsed as a float, and an error returned if parsing fails.
func NewProb(fstr string) (Prob, error) {
	if fstr == "*" {
		return MinProb, nil
	}

	f, err := strconv.ParseFloat(fstr, 64)
	if err != nil {
		return invalidProb,
			fmt.Errorf("Could not convert '%s' to a log probability: %s",
				fstr, err)
	}
	return Prob(f), nil
}

// Less returns true if `p1` represents a smaller probability than `p2`.
func (p1 Prob) Less(p2 Prob) bool {
	return p1 > p2
}

// IsMin returns true if the probability is minimal.
func (p Prob) IsMin() bool {
	return p == MinProb
}

// Ratio returns the log probability as a ratio in the range [0, 1].
// (This assumes that `p` is a log-odds score.)
func (p Prob) Ratio() float64 {
	if p.IsMin() {
		return 0.0
	}
	return math.Exp(-(float64(p)))
}

// Distance returns the distance between `p1` and `p2`.
func (p1 Prob) Distance(p2 Prob) float64 {
	return math.Abs(float64(p1) - float64(p2))
}

// String returns a string representation of the probability.
// When `p` is the minimum probability, then "*" is used.
// Otherwise, the full number is written.
func (p Prob) String() string {
	if p.IsMin() {
		return "*"
	}
	return fmt.Sprintf("%v", float64(p))
}

func (p Prob) MarshalJSON() ([]byte, error) {
	return json.Marshal(p.String())
}

func (p *Prob) UnmarshalJSON(bs []byte) error {
	var str string
	var err error
	if err := json.Unmarshal(bs, &str); err != nil {
		return err
	}
	if *p, err = NewProb(str); err != nil {
		return err
	}
	return nil
}

// HMMCat joins two HMMs together. The HMMs given are not modified.
// Both HMMs must have the same alphabet.
// The null emissions for the first HMM are used.
func HMMCat(h1, h2 *HMM) *HMM {
	nodes := make([]HMMNode, len(h1.Nodes)+len(h2.Nodes))
	copy(nodes, h1.Nodes)
	copy(nodes[len(h1.Nodes):], h2.Nodes)
	return &HMM{
		Nodes:    nodes,
		Alphabet: h1.Alphabet,
		Null:     h1.Null,
	}
}

// NewHMM creates a new HMM from a list of nodes, an ordered alphabet and a
// set of null probabilities (which may be nil).
func NewHMM(nodes []HMMNode, alphabet []Residue, null EProbs) *HMM {
	return &HMM{
		Nodes:    nodes,
		Alphabet: alphabet,
		Null:     null,
	}
}

// Slice returns a slice of the HMM given. A slice of an HMM returns only the
// HMM nodes (i.e., columns or match/delete states) in the region specified
// by the slice. Also, the transition probabilities of the last state are
// specially set: M->M = 0, M->I = *, M->D = *, I->M = 0, I->I = *, D->M = 0,
// and D->D = *.
// No other modifications are made.
func (hmm *HMM) Slice(start, end int) *HMM {
	nodes := make([]HMMNode, end-start)
	copy(nodes, hmm.Nodes[start:end])
	last := len(nodes) - 1
	nodes[last].Transitions.MM = 0
	nodes[last].Transitions.MI = MinProb
	nodes[last].Transitions.MD = MinProb
	nodes[last].Transitions.IM = 0
	nodes[last].Transitions.II = MinProb
	nodes[last].Transitions.DM = 0
	nodes[last].Transitions.DD = MinProb

	return &HMM{
		Nodes:    nodes,
		Alphabet: hmm.Alphabet,
		Null:     hmm.Null,
	}
}
