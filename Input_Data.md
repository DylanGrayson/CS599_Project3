# Phylogenetic/Metagenomic Binning
## Input Data Explanation

For Phylogenetic Binning there are two data sources that we need to pull from:
* A Sequence file
* A Reads file

### Sequence file

A sequence file contains a large amount of sequences...

A single sequence looks like this

```
>808338555-1423-Bacillus subtilis HJ5
ATCTTTTTCGTCTTTTTTTAGTATCCACAGAGGTTATCGACAAC ...
```

The length of some of these sequences is approximately 4 million characters long.

### Reads file


A single read looks like this

```
>SRR1029849.33351401 33351401 length=100
GAGCAGAATGTTAAACACATTAGCCCAGTCGCCGAAAACGGGGCCGATACTAGTCGACTGCATGCGAAATCTGTGGACACTGCCGGTGCGAATCGCGCGT
```

The length of a single read ranges from 100 to 250. 
