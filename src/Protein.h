/*
 * Protein.h
 *
 *  Created on: Sep 22, 2015
 *      Author: Chao
 */

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "Point.h"
#include <string>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
using namespace std;
class Protein {
public:
	Protein();
	Protein(string proteinName);
	string getProteinName();
	void setProteinName(string proteinName);
	int getRealSequenceLength();
	void setRealSequenceLength(int realSequenceLength);
	string getRealSequenceInfo() ;
	void setRealSequenceInfo(string realSeqeunceInfo);
	Point* getCAlpha_XYZ() ;
	void setCAlpha_XYZ(Point* CAlpha_XYZ);

	void freeProtein();
	virtual void loadProteinInfo(string)=0;

	void dispProteinInfo();

	static double calculateRMSD(Point*,Point*,int,int);


	virtual ~Protein();
protected:
	string proteinName;

	int realSequenceLength;
	string realSequenceInfo;
	Point* CAlpha_XYZ;




};

#endif /* PROTEIN_H_ */
