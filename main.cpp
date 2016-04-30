/*****
 * Dylan Grayson
 * Conner Swann
 * Brandon Paree
 *
 * CS599 Project 3
 *
 * Phylogenetic / Metagenomic Binning
******/


#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <algorithm>
#include <thread>
#include "Binning.hpp"

#define K 22
#define M 3

using namespace std;

map<string, Bucket*> getBucketList(char * filename);
void distributeReads(map<string, Bucket*> bList, char * seqFile, char* readFile);
bool extends(int readPos, string read, long int seqPos, ifstream* seqs);

//////////////////////////////////////////////////////
///////////////// MAIN ///////////////////////////////
int main(int argc, char* argv[]) {
	//Check that we have two inputs
	if (argc == 3) {
		//build our map of buckets from the sequences file
		map<string, Bucket*> b = getBucketList(argv[1]);

		//distribute our reads into the buckets
		distributeReads(b, argv[1], argv[2]);
	}else { //otherwise print help string
		cout << "USAGE: " << argv[0] << " <sequence file> <reads file>" << endl;
	}
	return 0;
}

void distributeReads(map<string, Bucket*> bList, char * seqFile, char* readFile) {
	ifstream seqs;
	seqs.open(seqFile);

	// for each bucket in our bucket list
	int counter = 0;

	// For threading we might be able to take the entire contents of this for loop,
	// and create a new function that takes in relevant data as inputs.
	// For each iteration of the thread we can start a new thread?

	for (map<string, Bucket*>::iterator bucket = bList.begin(); bucket != bList.end(); ++bucket) {

		//get the list of sequences from this bucket
		list<Sequence*> seqList = bucket->second->getSeqList();
		map<string, long int> *kmerList = new map<string, long int>;

		for (list<Sequence*>::iterator seq = seqList.begin(); seq != seqList.end(); ++seq) {
			long int location = (*seq)->getLocation();
			seqs.seekg(location);
			string line;

			getline(seqs, line);

			long int seqSize = line.size() - K + 1;

			//string *kmers = new string[seqSize];
			for (long int i = 0; i < seqSize; i++) {
				kmerList->insert(pair<string, long int>( line.substr(i, K), location+i) );
			}
			//kmerList->insert(kmerList->begin(), kmers, kmers + seqSize);
			//delete[] kmers;
		}
		//sort(kmerList->begin(), kmerList->end());
		string read;
		string label;
		long int readLocation;
		ifstream reads;

		reads.open(readFile);

		while ( getline(reads, read)) { //loop through every line of reads file
			if (read[0] != '>'){ //if this is the actual nucleotide sequence

				unsigned int sel = 0;
				string q = read.substr(0, K);

				while((sel+1)*K < read.size()) { //loop through every K sized chunk of the read
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
		reads.close();
		delete kmerList;

		cout << "Bucket " << counter << ") " << bucket->second->getReadCount() << endl;

		counter++;
	}
	seqs.close();
}

//extends* function, checks if a read extends w/ M mismatches
bool extends(int readPos, string read, long int seqPos, ifstream* seqs) {

	//get the "would be" position of the read in the sequence
	seqs->seekg(seqPos - readPos);

	//get our string read as a classic char * C string
	const char * cread = read.c_str();
	int readSize = read.size();
	char * segment = new char[readSize];

	//get the read-length segment of the corresponding sequence
	seqs->read(segment, readSize);
	int mm = 0;

	//for every position in the read/sequence
	for (int i = 0; i < readSize; i++) {
		if (i == readPos) {//if this is the K-sized segment from our seed match
			i += K; //skip ahead and don't check it again, OPTIMIZED!
		}else if (cread[i] != segment[i]) { //else if the current chars are mismatch
			mm++;	//add one mismatch
			if (mm > M) {	//if we blew our top on mismatches
				//get the hell outta there, we don't extend.
				delete[] segment;
				return false;
			}
		}
	}

	//we made it out unscathed, we extend.
	delete[] segment;
	return true;
}


//build map of speciesId to Bucket pointer from a fasta file
map<string, Bucket*> getBucketList(char * filename) {

	ifstream seqs;
	seqs.open(filename);
	map<string, Bucket*> bucketList;
	map<string, Bucket*>::iterator iter;
	string line;

	string sequenceId;
	string speciesId;
	string desc;
	int count;
	//Loop through the lines in the file, storing in line string.
	while ( getline(seqs, line) ){

		int i = 1;
		int dashCount = 0; //count dashes seen to know which part of id
		string nameComponent; //store strings between dashes

		if (line[0] == '>') { //if this is the label of the sequence
			while ( line[i] != ' ' && line[i] != '\0'){
				//while we aren't at the end of the label
				if (line[i] == '-') { //if we've reached the end of this component
					if (dashCount == 0){ //sequence id
						sequenceId = nameComponent;
					}else if (dashCount == 1) { //species id
						speciesId = nameComponent;
					}else if (dashCount == 2) { //description
						desc = nameComponent;
					}
					nameComponent = '\0'; //reset component
					dashCount++;
					i++;
				}
				nameComponent += line[i]; //build component
				i++;
			}
			iter = bucketList.find(speciesId);

			if (iter == bucketList.end()) { //if species id doesn't exists in map yet

				Bucket * b = new Bucket(speciesId); //make new bucket
				bucketList.insert(pair<string, Bucket*>(speciesId, b)); //and insert it in the map
			}
		}else{ //else this is the actual nucleotide sequence
				//The format of Fasta files guarantee that all variables
				//Will be set by the previous line loop.
			long int currLoc = seqs.tellg();
			//current location is at the end of the sequence, so we set
			//it to the beginning
			currLoc -= (line.size()+1)*sizeof(char);

			//create new sequence
			Sequence* s = new Sequence(sequenceId, speciesId, currLoc, desc);
			//and add it to the bucket
			bucketList.at(speciesId)->insertSequence(s);
		}
		count++;
	}
	seqs.close();

	return bucketList;
}
