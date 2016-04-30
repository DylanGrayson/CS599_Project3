# Phylogenetic/Metagenomic Binning
## Implementation Analysis

### Notes on Complexity

* Complexity is highly dependent on a number of factors:
  * Number of buckets
  * Size of data set


**For every bucket that exists, we need to go through this process**
```c++
while ( getline(reads, read)) {
  //if this is the actual nucleotide sequence
  if (read[0] != '>'){

    ...

    //loop through every K sized chunk of the read
    while((sel+1)*K < read.size()) {
      //try to find the current read in the kmer list
      map<string, long int>::iterator kmer = kmerList->find(q);

      //if it exists and it extends*
      if (kmer != kmerList->end() && extends(sel*K, read, kmer->second, &seqs)) {
        //we create the read object and insert it into the current bucket and break.
        Read * r = new Read(label.substr(0, label.find(" ")), readLocation);
        bucket->second->insertRead(r);
        break;
      }
      q = read.substr(++sel*K, K);
    }
  }else { // this is the sequence identifier
    //getting our location+1 in the file is the address of the read we store in the object
    readLocation = reads.tellg();
    readLocation += sizeof(char);
    label = read;
  }
}
```
