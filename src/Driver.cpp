//============================================================================
// Name        : LocalAlignmentTMScore.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstring>
#include <iostream>
#include <fstream>
#include "Protein.h"
#include "GroundTruthRoot.h"
#include "HitProtein.h"
using namespace std;
int getLength(string A) {
	int length = 0;
	for (int i = 0; i < A.length(); i++) {
		if (A[i] != '-') {
			length++;
		}
	}
	return length;
}

string getSeq(string A) {
	string seq;
	for (int i = 0; i < A.length(); i++) {
		if (A[i] == '-') {
			continue;
		}else{
			seq+=A[i];
		}
	}
	return seq;
}
string changeName(char name){
	if(name=='A'){
		return "ALA";
	}else if(name=='R'){
		return "ARG";
	}else if(name=='N'){
		return "ASN";
	}else if(name=='D'){
		return "ASP";
	}else if(name=='C'){
		return "CYS";
	}else if(name=='Q'){
		return "GLN";
	}else if(name=='E'){
		return "GLU";
	}else if(name=='G'){
		return "GLY";
	}else if(name=='H'){
		return "HIS";
	}else if(name=='I'){
		return "ILE";
	}else if(name=='L'){
		return "LEU";
	}else if(name=='K'){
		return "LYS";
	}else if(name=='M'){
		return "MET";
	}else if(name=='F'){
		return "PHE";
	}else if(name=='P'){
		return "PRO";
	}else if(name=='S'){
		return "SER";
	}else if(name=='T'){
		return "THR";
	}else if(name=='W'){
		return "TRP";
	}else if(name=='Y'){
		return "TYR";
	}else if(name=='V'){
		return "VAL";
	}else if(name=='B'){
		return "ASX";
	}else if(name=='Z'){
		return "GLX";
	}
}
int main(int argc, char* argv[]) {
string experimentLocation("/home/cf797/test/3DCOMBExperiment/HardCase/");
string inputFile(experimentLocation);
inputFile += argv[1];
inputFile += "/input.aln";
FILE* fptr = fopen((char*) inputFile.c_str(), "r");
if (fptr == NULL) {
	cout << "input file: " << inputFile << " can't open" << endl;
} else {
	int lineLength = 5000;
	char line[lineLength];
	string A;
	string B;

	while (fgets(line, lineLength, fptr) != NULL) {

		//cout << line << endl;
		char str[300];
		if (strstr(line, argv[1]) != NULL) {
			sscanf(line, "%*s %s", str);
			A += str;
		} else if (strstr(line, argv[2]) != NULL) {
			sscanf(line, "%*s %s", str);
			B += str;
		}
	}
	cout << A << endl;
	cout << endl;
	cout << B << endl;
	string Aseq=getSeq(A);
	string Bseq=getSeq(B);
	HitProtein protein(argv[2]);
	protein.loadProteinInfo("/home/lihongb/DATABASE/DBInfo/");

	int queryStart = 1;
	int queryEnd = getLength(A);
	string queryPart = A;
	int subjectStart = 1;
	int subjectEnd = getLength(B);
	string subjectPart = B;
	Point* fetchPart = protein.fetchSubjectAlignedPart3DPointsForQuery(
			queryStart, queryEnd, queryPart, subjectStart, subjectEnd,
			subjectPart);

	string outputFile(experimentLocation);
	outputFile+=argv[1];
	outputFile += "/";
	outputFile += argv[1];
	outputFile += "_";
	outputFile += argv[2];
	outputFile += ".pdb";
	cout<<outputFile<<endl;
	FILE *pFile=fopen((char*)outputFile.c_str(),"w");
	cout<<Aseq<<endl;
	for(int i=0;i<queryEnd-queryStart+1;i++){
		cout<<changeName(Aseq[i])<<","<<fetchPart[i].getX()<<","<<fetchPart[i].getY()<<","<<fetchPart[i].getZ()<<endl;
	}
	for (int i = 0; i < queryEnd - queryStart + 1; i++) {

		if(fetchPart[i].getX()==10000){
			continue;
		}
		string aminoAcid(changeName(Aseq[i]));
		fprintf(pFile,"ATOM  %5d  CA  %3s  %4d    %8.3f%8.3f%8.3f\n",i+1,(char*)aminoAcid.c_str(),i+1,fetchPart[i].getX(),fetchPart[i].getY(),fetchPart[i].getZ());

	}
	fprintf(pFile,"TER\n");
	fprintf(pFile,"END\n");
	fclose(pFile);
//protein.freeProtein();

	/*
	 GroundTruthRoot groundTruthRoot("T0821");
	 groundTruthRoot.setRealSequenceLength(255);
	 groundTruthRoot.loadProteinInfo("/home/cf797/test/DATABASE/PDB/");

	 cout<<"RMSD is:"<<Protein::calculateRMSD(fetchPart,groundTruthRoot.getCAlpha_XYZ(),queryStart,queryEnd);
	 */

	return 0;
}
}
