# MSAcleaner
a simple python script that takes in a FASTA formatted MSA and, given a set of 'reference' sequences, deletes any sequences that contribute to insertions given a fractional abundance

the output is then an un-aligned FASTA file ready for re-alignment

## example

```
python MSAcleaner.py -i input.msa -ref reference.seqIDs -fxn 0.01 -o output.aa
```

## to do:

I think the fractional abundance counter needs checking, ```-fxn 0``` seems to work, and low ```-fxn``` values are ok, but something goes wrong when there are more sequences? I think its to do with rounding errors
