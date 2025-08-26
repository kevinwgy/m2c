#include <math.h>
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {

	if(argc<2) {
		std::cout << "Usage: [binary] <surf-topo> <angle degs> <last-lmpc-id> <out-file>.\n";
		exit(-1);
	}

	double AngleOfSurface = atof(argv[2]);

	// read topo and find unique nodes
	int nNodesPerElement = 3; // assuming triangulated surface
	int *data = new int[nNodesPerElement];
	set<int> lmpcNodeSet;

	ifstream input(argv[1], ios::in);
	input.ignore(512, '\n');

	string line;
	while(getline(input, line)) {
		int eleId, eleType;
		istringstream is(line);
		is >> eleId >> eleType;
		if(eleType!=3) {
			std::cerr << "This utility expects a triangulated surface topology.\n";
			exit(-1);
		}

		for(int i=0; i<nNodesPerElement; ++i)
			is >> data[i];

		for(int i=0; i<nNodesPerElement; ++i)
			lmpcNodeSet.insert(data[i]);
	}

	std::cout << "Found " << lmpcNodeSet.size() << " unique nodes.\n";
	input.close();
	delete [] data;

	ofstream output(argv[4], ios::out);
	
	// write LMPC constraints of the form
	// 	sin(angle)*u_y + cos(angle)*u_z = 0.0
	output << "LMPC\n";
	int lmpcCounter = atoi(argv[3]);
	AngleOfSurface *= M_PI/180; // convert to radians
	double cosa = (AngleOfSurface < 0) ? -cos(AngleOfSurface) :  cos(AngleOfSurface);
	double sina = (AngleOfSurface < 0) ?  sin(AngleOfSurface) : -sin(AngleOfSurface);
	for(auto i=lmpcNodeSet.begin(); i!=lmpcNodeSet.end(); ++i) {
		output << setw(6) << ++lmpcCounter << setw(6) << 0 << "\n";
		output << setw(6) << *i << setw(6) << 2 << setw(16) << setprecision(12) << sina << "\n"; 
		output << setw(6) << *i << setw(6) << 3 << setw(16) << setprecision(12) << cosa << "\n"; 
	}	

	output.close();

	return 0;

}
