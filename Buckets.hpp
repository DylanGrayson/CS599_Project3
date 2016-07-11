/*
Class Definitions for:
 - Sequence
 - Read
 - Bucket
*/

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
		string getSequence() {
			return this->sequence;
		}
		Read(string rId, string seq) {
			readId = rId;
			sequence = seq;
		}
	private:
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
		int getReadCount() {
			return readList.size();
		}
		list<Sequence*> getSeqList() {
			return seqList;
		}
		void insertRead(Read* read) {
			readList.push_back(read);
		}
	private:
		list<Sequence*> seqList;
		list<Read*> readList;
};
