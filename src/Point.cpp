/*
 * Point.cpp
 *
 *  Created on: Sep 22, 2015
 *      Author: Chao
 */

#include "Point.h"

Point::Point() {
	// TODO Auto-generated ructor stub

}

Point::Point(double _x, double _y, double _z) {
	setX(_x);
	setY(_y);
	setZ(_z);
}

Point::~Point() {
	// TODO Auto-generated destructor stub
}


double Point::getX() {
	return x;
}

void Point::setX(double x) {
	this->x = x;
}

double Point::getY() {
	return y;
}

void Point::setY(double y) {
	this->y = y;
}

double Point::getZ() {
	return z;
}

void Point::setZ(double z) {
	this->z = z;
}
