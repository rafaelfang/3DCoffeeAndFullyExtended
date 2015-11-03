/*
 * HitProtein.cpp
 *
 *  Created on: Sep 26, 2015
 *      Author: Chao
 */

#include "HitProtein.h"



HitProtein::HitProtein():Protein() {
	// TODO Auto-generated constructor stub

}
HitProtein::HitProtein(string _proteinName) :
		Protein(_proteinName) {

}

Point* HitProtein::fetchSubjectAlignedPart3DPointsForQuery(int queryStart,
		int queryEnd, string queryPart, int subjectStart, int subjectEnd,
		string subjectPart) {

	int sizeOfAlignment = subjectPart.size();
	//cout << "align size" << sizeOfAlignment << endl;
	Point* subjectPoints = (Point*) malloc(sizeof(Point) * sizeOfAlignment);

	int numOfGapsInSubjectPart = 0;
	for (int i = 0; i < sizeOfAlignment; i++) {
		if (subjectPart[i] == '-') {
			Point bigNum(10000, 10000, 10000);
			subjectPoints[i] = bigNum;

			numOfGapsInSubjectPart++;
		} else {

			subjectPoints[i] = CAlpha_XYZ[subjectStart + i
					- numOfGapsInSubjectPart - 1];

		}

	}
	/*
	cout << "subjectPoints:" << endl;
	for (int i = 0; i < sizeOfAlignment; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << subjectPoints[i].getX() << " ";
	}
	cout << endl;
	for (int i = 0; i < sizeOfAlignment; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << subjectPoints[i].getY() << " ";
	}
	cout << endl;
	for (int i = 0; i < sizeOfAlignment; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << subjectPoints[i].getZ() << " ";
	}

	cout << endl;*/
	int numOfGapsInQueryPart = 0;
	for (int i = 0; i < queryPart.size(); i++) {
		if (queryPart[i] == '-') {
			numOfGapsInQueryPart++;
		}
	}

	int candidateQueryAlignedPartLength = queryEnd - queryStart + 1;
	Point* candidateQueryAlignedPart3DCoords = (Point*) malloc(
			sizeof(Point) * candidateQueryAlignedPartLength);

	int gapOffset = 0;
	for (int i = 0; i < queryPart.size(); i++) {
		if (queryPart[i] == '-') {
			gapOffset++;
		} else {
			candidateQueryAlignedPart3DCoords[i - gapOffset] = subjectPoints[i];
		}

	}
/*
	cout<<"query Part All XYZs:"<<endl;
	for (int i = 0; i < sizeOfAlignment - numOfGapsInQueryPart; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << candidateQueryAlignedPart3DCoords[i].getX() << " ";

	}
	cout << endl;

	for (int i = 0; i < sizeOfAlignment - numOfGapsInQueryPart; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << candidateQueryAlignedPart3DCoords[i].getY() << " ";

	}
	cout << endl;
	for (int i = 0; i < sizeOfAlignment - numOfGapsInQueryPart; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << candidateQueryAlignedPart3DCoords[i].getZ() << " ";

	}

	cout << endl;
*/
	return candidateQueryAlignedPart3DCoords;

}
void HitProtein::loadProteinInfo(string proteinLocation) {
	//open protein db file, read and store protein infomation
	string proteinFile(proteinLocation);
	proteinFile += proteinName;
	proteinFile += ".db";
	FILE* fptr = fopen((char*) proteinFile.c_str(), "r");
	if (fptr == NULL) {
		cout << "input file: " << proteinFile << " can't open" << endl;
	} else {
		int lineLength = 5000;
		char line[lineLength];
		while (fgets(line, lineLength, fptr) != NULL) {
			if ((strstr(line, ">Real Sequence Info:") != NULL)) {

				fgets(line, lineLength, fptr);

				int numberOfChars = 0;
				while (line[numberOfChars] != '\0') {

					numberOfChars++;
				}
				numberOfChars--;

				string _realSequenceInfo(line);

				setRealSequenceInfo(_realSequenceInfo);
				//cout << realSequenceInfo << endl;
				setRealSequenceLength(numberOfChars);
				//cout << realSequenceLength << endl;
			}
			if ((strstr(line, ">Ca XYZ:") != NULL)) {

				double Xs[realSequenceLength];
				double Ys[realSequenceLength];
				double Zs[realSequenceLength];

				double temp;
				for (int i = 0; i < realSequenceLength; i++) {
					fscanf(fptr, "%lf", &temp);
					Xs[i] = temp;
					//if (i % 10 == 0) {
					// cout << endl;
					// }
					// cout << Xs[i] << " ";
				}
				//cout << endl;

				for (int i = 0; i < realSequenceLength; i++) {
					fscanf(fptr, "%lf", &temp);
					Ys[i] = temp;
					//if (i % 10 == 0) {
					// cout << endl;
					// }
					// cout << Ys[i] << " ";
				}
				//cout << endl;
				for (int i = 0; i < realSequenceLength; i++) {
					fscanf(fptr, "%lf", &temp);
					Zs[i] = temp;
					//if (i % 10 == 0) {
					// cout << endl;
					// }
					// cout << Zs[i] << " ";
				}
				//cout << endl;
				CAlpha_XYZ = (Point*) malloc(
						sizeof(Point) * realSequenceLength);

				for (int i = 0; i < realSequenceLength; i++) {
					Point point(Xs[i], Ys[i], Zs[i]);
					CAlpha_XYZ[i] = point;
				}

			}
		}
	}

}
HitProtein::~HitProtein() {
	// TODO Auto-generated destructor stub
}


