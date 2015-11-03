/*
 * HitProtein.h
 *
 *  Created on: Sep 26, 2015
 *      Author: Chao
 */

#ifndef HITPROTEIN_H_
#define HITPROTEIN_H_

#include "Protein.h"

class HitProtein: public Protein {
public:
	HitProtein();
	HitProtein(string);
	virtual ~HitProtein();
	Point* fetchSubjectAlignedPart3DPointsForQuery(int, int, string, int, int,
			string);
	void loadProteinInfo(string);
};

#endif /* HITPROTEIN_H_ */
