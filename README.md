# CS599_Project3

##Project #3

###Phylogenetic / Metagenomic Binning

####Project Intro:

Recent advances in High Throughput Sequencing (HTS) technologies have made shotgun sequencing of complex environmental samples economically feasible. In theory, such approaches are capable of characterizing the entirety of a sample’s composition, simultaneously identifying all present eukaryotes (plants, animals, and fungi), prokaryotes, and viruses. The first and foremost challenge of metagenomic characterization of unknown samples is that of ‘binning’ sequence fragments (reads) to associate them with species potentially present in the sample. For this project, your goal will be to identify which sequence reads can be associated with which organisms with up to a predefined number of tolerated mismatches.

At its core, this is a counting exercise – how many sequence read fragments are suggesting presence of any given candidate organism. Any given read can simultaneously belong to multiple species bins or to none (with up to M mismatches).

####Inputs:

* High Throughput Sequencing (HTS) reads from a real metagenomic re-sequencing experiments. You will be required to test your pipeline on multiple read sets (provided by instructor) ranging in size from 1 million 50-mers to 100 million 50 mers.
* A set of bins (k = 1,…,m), each corresponding to a bacterial species, with multiple reference genomes in each.
* Total number of tolerated mismatches.

####Goals:
* For all reads/queries find the first (if it exists) alignment between the read and any sequence S in each bin with up to M mismatches.  

* You must return the total number of reads belonging to each bin (distinct species) and the list of bins aligning to each species. Alignment quality does not matter as long as the alignment is with up to M mismatches.  This means that if M = 3, a perfect alignment is not any ‘better’ than the one with 2 mismatches.

####Additional requirements:

* Unless explicitly approved by the instructor, the core the pipeline must be written in C or C++.  Python wrappers to handle data input and output are OK.

* You must be able to ‘escape’ the alignment process once the first alignment between a read and a bin is achieved – there is no point in looking for additional alignments (hence binning).

* You must provide extensive documentation on all classes / functions that you create, as well as documentation for all external libraries you are incorporating into your project.

* You must provide benchmarks and time /space estimates on the performance of your pipeline.

#### Compilation Instructions

```
g++ -o <output file> main.cpp -std=c++11 -lgomp -pthread
```

#### Usage Instructions

```
./<output file> <sequence file> <reads file>
```
