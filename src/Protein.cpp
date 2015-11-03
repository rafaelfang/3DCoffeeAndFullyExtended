/*
 * Protein.cpp
 *
 *  Created on: Sep 22, 2015
 *      Author: Chao
 */

#include "Protein.h"

Protein::Protein() {

}

Protein::Protein(string _proteinName) {
	setProteinName(_proteinName);

}

string Protein::getProteinName() {
	return proteinName;
}

void Protein::setProteinName(string proteinName) {
	this->proteinName = proteinName;
}

int Protein::getRealSequenceLength() {
	return realSequenceLength;
}

void Protein::setRealSequenceLength(int realSequenceLength) {
	this->realSequenceLength = realSequenceLength;
}

string Protein::getRealSequenceInfo() {
	return realSequenceInfo;
}

void Protein::setRealSequenceInfo(string realSequenceInfo) {
	this->realSequenceInfo = realSequenceInfo;
}

Point* Protein::getCAlpha_XYZ() {
	return CAlpha_XYZ;
}

void Protein::setCAlpha_XYZ(Point* _CAlpha_XYZ) {
	CAlpha_XYZ = _CAlpha_XYZ;
}

double Protein::calculateRMSD(Point* fetchSubjectAlignedPart3DPointsForQuery,Point*GroundTruthRootCAlpha_XYZ,int queryStart,int queryEnd){
	int n=0;
	double temp=0;
	for(int i=0;i<queryEnd-queryStart+1;i++){
		if(fetchSubjectAlignedPart3DPointsForQuery[i].getX()==10000){
			continue;
		}else{

			double diffX=fetchSubjectAlignedPart3DPointsForQuery[i].getX()-GroundTruthRootCAlpha_XYZ[queryStart+i-1].getX();
			//cout<<"HItProtein X: "<<fetchSubjectAlignedPart3DPointsForQuery[i].getX()<<" GroundX: "<<GroundTruthRootCAlpha_XYZ[queryStart+i-1].getX()<<endl;
			double diffY=fetchSubjectAlignedPart3DPointsForQuery[i].getY()-GroundTruthRootCAlpha_XYZ[queryStart+i-1].getY();
			double diffZ=fetchSubjectAlignedPart3DPointsForQuery[i].getZ()-GroundTruthRootCAlpha_XYZ[queryStart+i-1].getZ();
			temp+=(diffX*diffX+diffY*diffY+diffZ*diffZ);
			n++;
		}
	}
//cout<<"temp: "<<temp<<"n "<<n<<endl;
	temp=temp/n;
	temp=sqrt(temp);


	return temp;
}



void Protein::dispProteinInfo() {
	cout << "Protein " << proteinName << endl;
	cout << "The real sequence length is " << realSequenceLength << endl;
	cout << "The real seqeunce information is " << endl;
	cout << realSequenceInfo << endl;
	cout << "The CAlpha_XYZ is as followed:" << endl;

	for (int i = 0; i < realSequenceLength; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << CAlpha_XYZ[i].getX() << " ";

	}

	cout << endl;

	for (int i = 0; i < realSequenceLength; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << CAlpha_XYZ[i].getY() << " ";

	}

	cout << endl;

	for (int i = 0; i < realSequenceLength; i++) {
		if (i % 10 == 0) {
			cout << endl;
		}
		cout << CAlpha_XYZ[i].getZ() << " ";

	}

	cout << endl;
}

void Protein::freeProtein() {
	free(CAlpha_XYZ);
}
Protein::~Protein() {
// TODO Auto-generated destructor stub
}

