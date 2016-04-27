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
#include <vector>
#include <algorithm>


#define K 12
#define M 3

using namespace std;


//Sequence Class contains a long int sequenceLoc that gets passed to 
//file.seekg() to get the actual sequence.
class Sequence {
	public:
		string sequenceId;
		string	speciesId;
		string description;
		long int getLocation() {
			return this->sequenceLoc;
		}
		Sequence(string seqId, string species, long int location, string desc) {
			sequenceId = seqId;
			speciesId = species;
			sequenceLoc = location;
			description = desc;
		}
	private:
		long int sequenceLoc;
};

//Read class
class Read {
	public:
		string readId;
		long int getLocation() {
			return this->readLoc;
		}
		Read(string rId, long int location) {
			readId = rId;
			readLoc = location;
		}
	private:
		long int readLoc;
};
		
//Bucket class contains a sequence list and a reads list
class Bucket {
	public:
		string speciesId;
		Bucket(string species) {
			speciesId = species;
		}
		int insertSequence(Sequence * seq) {
			seqList.push_back(seq);
			return seqList.size();
		}
		int getSeqCount() {
			return seqList.size();
		}
		list<Sequence*> getSeqList() {
			return seqList;
		}
	private:
		list<Sequence*> seqList;
		list<Read*> readList;
};

map<string, Bucket*> getBucketList(char * filename);
void distributeReads(map<string, Bucket*> bList, char * seqFile, char* readFile);
bool binarySearch(vector<string> *arr, vector<string>::iterator start, vector<string>::iterator end, string query);

		

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
	for (map<string, Bucket*>::iterator bucket = bList.begin(); bucket != bList.end(); ++bucket) {
		//get the list of sequences from this bucket
		list<Sequence*> seqList = bucket->second->getSeqList();
		vector<string> *kmerList = new vector<string>;
		for (list<Sequence*>::iterator seq = seqList.begin(); seq != seqList.end(); ++seq) {
			seqs.seekg((*seq)->getLocation());
			string line;
			getline(seqs, line);
			long int seqSize = line.size() - K + 1;
			string *kmers = new string[seqSize];
			for (long int i = 0; i < seqSize; i++) {
				kmers[i] = line.substr(i, K);
			}
			kmerList->insert(kmerList->begin(), kmers, kmers + seqSize);
			delete[] kmers;
		}
		sort(kmerList->begin(), kmerList->end());
		string read;
		string label;
		ifstream reads;
		reads.open(readFile);
		while ( getline(reads, read)) {
			if (read[0] != '>'){
				unsigned int sel = 0;
				string q = read.substr(read.size() - K, K);
				while((sel+1)*K < read.size()) {
					if(binarySearch(kmerList, kmerList->begin(), kmerList->end(), q)) {
						//cout << q << "-" << sel <<  " " << label << endl;
						break;
					}
					q = read.substr(sel*K, K);
					sel++;
				}
			}else {
				label = read;
			}
		}
		reads.close();
		delete kmerList;
		cout << counter << ") " << bucket->second->getSeqCount() << endl;
		counter++;
	}
	seqs.close();
}

bool binarySearch(vector<string> *arr, vector<string>::iterator start, vector<string>::iterator end, string query) {
	int length = end - start;
	if (length <= 1) {
		return (arr->at(0) == query);
	}else {
		int half = length / 2;
		string val = arr->at(half);
		if (query == val) {
			return true;
		}else if (query < val) {
			return binarySearch(arr, start, start+half, query);
		}else {
			return binarySearch(arr, start+half, end, query);
		}
	}
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
	
