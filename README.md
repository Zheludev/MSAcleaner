# MSAcleaner
a simple python script that takes in a FASTA formatted MSA and, given a set of 'reference' sequences, deletes any sequences that contribute to insertions given a fractional abundance

## example

```
python MSAcleaner.py -i input.msa -ref reference.seqIDs -fxn 0.01 -o output.aa
```
