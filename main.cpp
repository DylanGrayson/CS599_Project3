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
#include <stdlib.h>
using namespace std;

class Sequence {
	public:
		int readId;
		int	speciesId;
		string description;
		
		unsigned int getLocation() {
			return this->sequenceLoc;
		}
	private:
		unsigned int sequenceLoc;
};

class Read {
	public:
		int readId;
		string sequence;
};
		

class Bucket {
	public:
		int speciesId;
	private:
		Read * readList[];
		Sequence * seqList[];
};
		

int main(int argc, char* argv[]) {
	//All this should be deleted, just messing around with file I/O
	if (argc == 3) {
		
		ifstream seqs;
		seqs.open(argv[1]);
		string line;
		string names[71099][3];
		int count;
		while ( getline(seqs, line) ){
			int i = 1;
			int comp_count = 0;
			string name_component;
			if (line[0] == '>') {
				while ( line[i] != ' ' ){
					if (line[i] == '-') {
						names[count][comp_count++] = name_component;
						name_component = '\0';
					}
					name_component += line[i];
					i++;
				}
				cout << names[count][2] << endl;
			}
			count++;
		}
		seqs.close();
	}
	return 0;
}
