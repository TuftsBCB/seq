/*
Package seq provides common types and operations for dealing with biological
sequence data, with a bias toward amino acid sequences. Types includes
sequences, profiles, multiple sequence alignments and HMMs. Operations include
sequence alignment (currently only Needleman-Wunsch global alignment), building
frequency profiles with background probabilities and an implementation of the
Viterbi algorithm to find the probability of the most likely alignment of a
sequence to an HMM.

This package is currently a "kitchen sink" of operations on biological
sequences. It isn't yet clear (to me) whether it should remain a kitchen sink.
*/
package seq
