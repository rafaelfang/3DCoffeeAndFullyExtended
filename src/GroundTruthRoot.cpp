/*
 * GroundTruthRoot.cpp
 *
 *  Created on: Sep 27, 2015
 *      Author: Chao
 */

#include "GroundTruthRoot.h"

GroundTruthRoot::GroundTruthRoot() {
	// TODO Auto-generated constructor stub

}
GroundTruthRoot::GroundTruthRoot(string _proteinName) :
		Protein(_proteinName) {

}

void GroundTruthRoot::loadProteinInfo(string proteinLocation) {
	string groundTruthRootPDB(proteinLocation);
	groundTruthRootPDB += proteinName;
	groundTruthRootPDB += ".pdb";
	FILE* fptr = fopen((char*) groundTruthRootPDB.c_str(), "r");
	if (fptr == NULL) {
		cout << "input file: " << groundTruthRootPDB << " can't open" << endl;
	} else {
		int lineLength = 5000;
		char line[lineLength];
		CAlpha_XYZ=(Point*) malloc(sizeof(Point)*realSequenceLength);
		int index=0;
		while (fgets(line, lineLength, fptr) != NULL) {
			if ((strstr(line, "CA") != NULL)) {
				char info0[4];
				int info1;
				char info2[3], info3[4];
				int info4;
				double info5,info6,info7,info8,info9;
				char info10[1],info11[2];
				sscanf(line,"%s %d %s %s %d %lf %lf %lf %lf %lf %s %s",
						info0, &info1, info2, info3, &info4, &info5,
						&info6, &info7, &info8, &info9, info10, info11);
				Point p(info5, info6, info7);
				CAlpha_XYZ[index]=p;
				index++;

			}
		}
	}
}
GroundTruthRoot::~GroundTruthRoot() {
	// TODO Auto-generated destructor stub
}

