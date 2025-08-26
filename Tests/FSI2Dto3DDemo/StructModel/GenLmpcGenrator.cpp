#include <math.h>
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

void getNodalCoordinates(istream &input, vector<vector<double> > &coords);
void getSurfaceElements(istream &input, vector<vector<int> > &elem);
void getUniqueNodesAndNormals(const vector<vector<double> > &coords, const vector<vector<int> > &elems, 
			      set<int> &uniqnodes, vector<vector<double> > &normal);
void writeLmpcInclude(ostream &output, const set<int> &uniqnodes, const vector<vector<double> > &normal);	

int main(int argc, char* argv[]) {

	if(argc<3) {
		std::cout << "Usage: [binary] <msh-file> <surftopo-file> <out-file>.\n";
		exit(-1);
	}

	ifstream meshInput(argv[1], ios::in);
	ifstream surfInput(argv[2], ios::in);
	ofstream output(argv[3], ios::in);

	vector<vector<double> > nodeCoordinates;
	vector<vector<int> > surfElements;
	vector<vector<double> > lmpcNodeNormals;
	set<int> lmpcNodes;	

	getNodalCoordinates(meshInput, nodeCoordinates);
	getSurfaceElements(surfInput, surfElements);

	getUniqueNodesAndNormals(nodeCoordinates, surfElements, lmpcNodes, lmpcNodeNormals);
	//writeLmpcInclude(output, lmpcNodes, lmpcNodeNormals);	

	meshInput.close();
	surfInput.close();
	output.close();
	return 0;
}

void
getNodalCoordinates(istream &input, vector<vector<double> > &coords)
{

	int initialSize = 1000;
	if(coords.empty())
		coords.assign(initialSize, vector<double>(3, 0.0));
	
	input.ignore(512, '\n');
	int numNodes;
	double *data = new double[3];
	for(numNodes=0; numNodes<INT_MAX; ++numNodes) {
		int nodeId;
		bool done=false;
		string line;
		getline(input, line);

		if(line.compare("TOPOLOGY") == 0) break;

		istringstream is(line);
		is >> nodeId;
		for(int i=0; i<3; ++i) {
			is >> data[i];
			if(is.fail()) {
				done = true;
				break;
			}
		}

		if(done) break;
		if(numNodes>=coords.size())
			coords.resize(coords.size()+initialSize, vector<double>(3, 0.0));
		for(int i=0; i<3; ++i)
			coords[numNodes][i] = data[i];

	}

	coords.resize(numNodes);
	delete [] data;
}

void 
getSurfaceElements(istream &input, vector<vector<int> > &elem)
{

	int initialSize = 1000;
	if(elem.empty())
		elem.assign(initialSize, vector<int>(3, -1));

	input.ignore(512, '\n');
	
	int numElems;
	int *data = new int[3];
	for(numElems=0; numElems<INT_MAX; ++numElems) {
		int eleId, eleType;
		bool done = false;
		string line;
		getline(input, line);
		istringstream is(line);
		is >> eleId >> eleType;

		if(eleType != 3) {
			cerr << "***Error: Utility expectes a triangulated surface.\n";
			exit(-1);
		}
			
		for(int i=0; i<3; ++i) {
			is >> data[i];
			if(is.fail()) {
				done = true;
				break;
			}
		}
		if(done) break;
		if(numElems>=elem.size())
			elem.resize(initialSize+elem.size(), vector<int>(3,-1));
		for(int i=0; i<3; ++i)
			elem[numElems][i] = data[i];

	}

	elem.resize(numElems);
	delete [] data;
}

void 
getUniqueNodesAndNormals(const vector<vector<double> > &coords, const vector<vector<int> > &elems, 
			 set<int> &uniqnodes, vector<vector<double> > &normal)
{

	// compute normals for all surface elements
	vector<vector<double> > elementNormals(elems.size(), vector<double>(3,0.0));
	for(int iele=0; iele<elems.size(); ++iele) {
		vector<double> A = coords[elems[iele][0]-1]; 	
		vector<double> B = coords[elems[iele][1]-1]; 	
		vector<double> C = coords[elems[iele][2]-1];
		
		vector<double> BA(3, 0.0), CA(3, 0.0);
		for(int i=0; i<3; ++i) {
			BA[i] = B[i] - A[i];
			CA[i] = C[i] - A[i];
		}
	
		// compute cross product
		elementNormals[iele][0] = BA[1]*CA[2]-BA[2]*CA[1];
		elementNormals[iele][1] = BA[2]*CA[0]-BA[0]*CA[2];
		elementNormals[iele][2] = BA[0]*CA[1]-BA[1]*CA[0];

		// normalize
		double norm = sqrt(elementNormals[iele][0]*elementNormals[iele][0]+
				   elementNormals[iele][1]*elementNormals[iele][1]+
				   elementNormals[iele][2]*elementNormals[iele][2]);
		elementNormals[iele][0] /= norm;
		elementNormals[iele][1] /= norm;
		elementNormals[iele][2] /= norm;
	}

	// find all unique nodes ids on the surface topology
	for(int iele=0; iele<elems.size(); ++iele) 
		for(int inode=0; inode<elems[iele].size(); ++inode) 
			uniqnodes.insert(elems[iele][inode]);

	cout << "Found: " << uniqnodes.size() << " unique nodes.\n";
	normal.assign(uniqnodes.size(), vector<double>(3, 0.0));

	int indexLocation=0;
	for(auto inode=uniqnodes.begin(); inode!=uniqnodes.end() ; ++inode, ++indexLocation) {
		// get all elements which are attached to the node
		int currentNode = *inode;
		set<int> attachedElements;
		// brute force approach
		for(int iele=0; iele<elems.size(); ++iele) {
			for(int i=0; i<elems[iele].size(); ++i) {
				if(currentNode == elems[iele][i]) {
					attachedElements.insert(iele);
					break;
				}
			}
		}

		if(attachedElements.size() == 0) {
			cerr << "***Error: Node " << currentNode << " does not belong to any surface element.\n";
			exit(-1);
		}

		for(auto iter=attachedElements.begin(); iter!=attachedElements.end(); ++iter) {
			normal[indexLocation][0] += elementNormals[*iter][0]/attachedElements.size();
			normal[indexLocation][1] += elementNormals[*iter][1]/attachedElements.size();
			normal[indexLocation][2] += elementNormals[*iter][2]/attachedElements.size();
		}

		// clean-up normals
		for(int j=0; j<3; ++j)
			if(abs(normal[indexLocation][j]) < 1e-3) normal[indexLocation][j] = 0.0;
	}
}

void 
writeLmpcInclude(ostream &output, const set<int> &uniqnodes, const vector<vector<double> > &normal)
{

}



















































