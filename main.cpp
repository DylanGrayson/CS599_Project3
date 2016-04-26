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
	private:
		list<Sequence*> seqList;
		list<Read*> readList;
};

map<string, Bucket*> getBucketList(char * filename);
		

int main(int argc, char* argv[]) {
	//All this should be deleted, just messing around with file I/O
	if (argc == 2) {
		map<string, Bucket*> b = getBucketList(argv[1]);
		for (map<string, Bucket*>::iterator iter = b.begin(); iter != b.end(); ++iter) {
			int seqCount = iter->second->getSeqCount();
			if (seqCount > 100) {
				cout << iter->first << " => " << seqCount << endl;
			}
		}
	}else {
		cout << "must pass in sequences file" << endl;
	}
	return 0;
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
			Sequence* s = new Sequence(sequenceId, speciesId, seqs.tellg(), desc);
			bucketList.at(speciesId)->insertSequence(s);
		}
		count++;
	}
	seqs.close();
	return bucketList;
}
	
