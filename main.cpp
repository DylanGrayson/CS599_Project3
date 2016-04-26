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

#define K 22

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
		string sequence;
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
	ifstream reads;
	reads.open(readFile);
	// for each bucket in our bucket list
	for (map<string, Bucket*>::iterator bucket = bList.begin(); bucket != bList.end(); ++bucket) {
		//get the list of sequences from this bucket
		list<Sequence*> seqList = bucket->second->getSeqList();
		string sequence;
		for (list<Sequence*>::iterator seq = seqList.begin(); seq != seqList.end(); ++seq) {
			seqs.seekg((*seq)->getLocation());
			string line;
			getline(seqs, line);
			sequence += line;
		}
		int seqSize = sequence.size();
		string* kmerList = new string[seqSize - K + 1];
		for (int i = 0; i < seqSize - K + 1; i++) {
			kmerList[i] = sequence.substr(i, K);
		}
		delete[] kmerList;
		cout << bucket->second->speciesId << endl;
	}
	seqs.close();
	reads.close();
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
	
