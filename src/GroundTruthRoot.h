/*
 * GroundTruthRoot.h
 *
 *  Created on: Sep 27, 2015
 *      Author: Chao
 */

#ifndef GROUNDTRUTHROOT_H_
#define GROUNDTRUTHROOT_H_

#include "Protein.h"

class GroundTruthRoot: public Protein {
public:
	GroundTruthRoot();
	GroundTruthRoot(string);
	void loadProteinInfo(string);
	virtual ~GroundTruthRoot();
};

#endif /* GROUNDTRUTHROOT_H_ */
