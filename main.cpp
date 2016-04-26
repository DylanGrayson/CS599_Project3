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

class Read {
	public:
		string readId;
		string sequence;
};
		

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
	//All this should be deleted, just messing around with file I/O
	if (argc == 3) {
		map<string, Bucket*> b = getBucketList(argv[1]);
		distributeReads(b, argv[1], argv[2]);
	}else {
		cout << "Usage: " << argv[0] << " <sequence file> <reads file>" << endl;
	}
	return 0;
}

void distributeReads(map<string, Bucket*> bList, char * seqFile, char* readFile) {
	ifstream seqs;
	seqs.open(seqFile);
	ifstream reads;
	reads.open(readFile);
	int poopoopoo = 0;
	for (map<string, Bucket*>::iterator bucket = bList.begin(); bucket != bList.end(); ++bucket) {
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
	while ( getline(seqs, line) ){
		int i = 1;
		int dashCount = 0;
		string nameComponent;
		if (line[0] == '>') {
			while ( line[i] != ' ' && line[i] != '\0'){
				if (line[i] == '-') {
					if (dashCount == 0){
						sequenceId = nameComponent;
					}else if (dashCount == 1) {
						speciesId = nameComponent;
					}else if (dashCount == 2) {
						desc = nameComponent;
					}
					nameComponent = '\0';
					dashCount++;
					i++;
				}
				nameComponent += line[i];
				i++;
			}
			iter = bucketList.find(speciesId);
			if (iter == bucketList.end()) {
				Bucket * b = new Bucket(speciesId);
				bucketList.insert(pair<string, Bucket*>(speciesId, b));
			}
		}else{
			long int currLoc = seqs.tellg();
			currLoc -= (line.size()+1)*sizeof(char);
			Sequence* s = new Sequence(sequenceId, speciesId, currLoc, desc);
			bucketList.at(speciesId)->insertSequence(s);
		}
		count++;
	}
	seqs.close();
	return bucketList;
}
	
