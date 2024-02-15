/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.25 $                                                              
// $Date: 2009-01-29 00:42:03 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MCTM.cpp,v $                                                                
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for MCTM.
//
// What: "@(#) MCTM.C, revA"

#include <string.h>

#include <MCTM.h>
#include <UniaxialMaterial.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>
#include <MaterialResponse.h>

#define sqrt2 1.41421356237310

void*
OPS_MCTM(void)
{
	NDMaterial* theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs != 10 && numArgs != 11) {
		opserr << "Want: nDMaterial MCTM $tag $m1 ~ $m9 <$m10>" << endln;
		return 0;
	}

	int tag;
	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
		return 0;
	}

	numArgs = OPS_GetNumRemainingInputArgs();
	int* iData = new int[numArgs];
		
//		new int[numArgs];
	if (OPS_GetInt(&numArgs, iData) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
		return 0;
	}

	UniaxialMaterial** concrete = new UniaxialMaterial* [9];
	UniaxialMaterial** steel = nullptr;

	for (int i = 0; i < 9; i++) {
		concrete[i] = nullptr;
		concrete[i] = OPS_getUniaxialMaterial(iData[i]);
	}

	for (int i = 0; i < 9; i++) {
		if (concrete[i] == nullptr)
		{
			opserr << "WARNING invalid integer tag: " << iData[i] << "UniaxialMaterial when building MCTM \n";
			return theMaterial;
		}
	}

	if (numArgs == 9)
		theMaterial = new MCTM(tag, concrete);
	if (numArgs == 10) {
		steel = new UniaxialMaterial* [1];
		steel[0] = OPS_getUniaxialMaterial(iData[9]);
		if (steel[0] == nullptr)
		{
			opserr << "WARNING invalid integer tag: " << iData[9] << "UniaxialMaterial when building MCTM \n";
			return theMaterial;
		}
		theMaterial = new MCTM(tag, concrete, steel);
		delete[] steel;
	}

	delete[] concrete;
	delete[] iData;

	return theMaterial;
}

Matrix MCTM::T(2, 2);
Matrix MCTM::Ttran(2, 2);

MCTM::MCTM(int tag, UniaxialMaterial** c)
	: NDMaterial(tag, MAT_TAG_MCTM), rho(0.0), steel(nullptr),
	eps(6), sig(6), K(6, 6), Ki(0), etruss(9), struss(9)
{
	concrete = new UniaxialMaterial* [9];

	for (int i = 0; i < 9; i++)
		concrete[i] = c[i]->getCopy();

	T.Zero();

	T(0, 0) = 1 / sqrt2;
	T(1, 0) = 1 / sqrt2;
	T(0, 1) = -1 / sqrt2;
	T(1, 1) = 1 / sqrt2;

	Ttran.addMatrixTranspose(0, T, 1);
}

MCTM::MCTM(int tag, UniaxialMaterial** c, UniaxialMaterial** s)
	: NDMaterial(tag, MAT_TAG_MCTM), rho(0.0),
	eps(6), sig(6), K(6, 6), Ki(0), etruss(10), struss(10)
{
	concrete = new UniaxialMaterial * [9];

	for (int i = 0; i < 9; i++)
		concrete[i] = c[i]->getCopy();

	steel = new UniaxialMaterial * [1];
	steel[0] = s[0]->getCopy();

	T.Zero();

	T(0, 0) = 1 / sqrt2;
	T(1, 0) = 1 / sqrt2;
	T(0, 1) = -1 / sqrt2;
	T(1, 1) = 1 / sqrt2;

	Ttran.addMatrixTranspose(0, T, 1);
}

MCTM::~MCTM()
{
	if (Ki != 0)
		delete Ki;
	for (int i = 0; i < 9; i++)
		delete concrete[i];
	delete[] concrete;

	if (steel != nullptr)
	{
		delete steel[0];
		delete[] steel;
	}
}

double
MCTM::getRho()
{
	return rho;
}

NDMaterial*
MCTM::getCopy(const char* type)
{
	NDMaterial* theMat = new MCTM(this->getTag(), concrete, steel);

	return theMat;
}

int
MCTM::setTrialStrain(const Vector& v)
{
	eps = v; // xx,yy,zz,xy,yz,xz
	K.Zero();
	sig.Zero();

	concrete[0]->setTrialStrain(v(0));
	sig(0) = concrete[0]->getStress();
	K(0, 0) = concrete[0]->getTangent();
	etruss(0) = v(0);
	struss(0) = sig(0);

	concrete[1]->setTrialStrain(v(1));
	sig(1) = concrete[1]->getStress();
	K(1, 1) = concrete[1]->getTangent();
	etruss(1) = v(1);
	struss(1) = sig(1);

	concrete[2]->setTrialStrain(v(2));
	sig(2) = concrete[2]->getStress();
	K(2, 2) = concrete[2]->getTangent();
	etruss(2) = v(2);
	struss(2) = sig(2);
	if (steel != nullptr)
	{
		steel[0]->setTrialStrain(v(0));
		sig(2) += steel[0]->getStress();
		K(2, 2) += steel[0]->getTangent();
		etruss(9) = v(2);
		struss(9) = sig(2);
	}

	Matrix epsg(2, 2), epsl(2, 2), sigg(2, 2), sigl(2, 2), tgtg(2, 2), tgtl(2, 2);

	epsl.Zero();
	sigl.Zero();
	tgtl.Zero();
	epsg(0, 0) = v(0) - (v(0) + v(1)) / 2;
	epsg(1, 1) = v(1) - (v(0) + v(1)) / 2;
	epsg(0, 1) = v(3);
	epsg(1, 0) = v(3);
	epsl = Ttran * epsg * T;
	etruss(3) = epsl(0, 0);
	etruss(4) = epsl(1, 1);
	concrete[3]->setTrialStrain(epsl(0,0));
	concrete[4]->setTrialStrain(epsl(1,1));
	struss(3) = concrete[3]->getStress();
	struss(4) = concrete[4]->getStress();
	sigl(0, 0) = struss(3);
	sigl(1, 1) = struss(4);
	tgtl(0, 0) = concrete[3]->getTangent();
	tgtl(1, 1) = -concrete[4]->getTangent();
	sigg = T * sigl * Ttran;
	tgtg = T * tgtl * Ttran;
	sig(0) += sigg(0, 0);
	sig(1) += sigg(1, 1);
	sig(3) += sigg(0, 1);
	K(0, 0) += tgtg(0, 0);
	K(1, 1) += tgtg(1, 1);
	K(3, 3) += tgtg(0, 1);

	epsl.Zero();
	sigl.Zero();
	tgtl.Zero();
	epsg(0, 0) = v(1) - (v(1) + v(2)) / 2;
	epsg(1, 1) = v(2) - (v(1) + v(2)) / 2;
	epsg(0, 1) = v(4);
	epsg(1, 0) = v(4);
	epsl = Ttran * epsg * T;
	etruss(5) = epsl(0, 0);
	etruss(6) = epsl(1, 1);
	concrete[5]->setTrialStrain(epsl(0, 0));
	concrete[6]->setTrialStrain(epsl(1, 1));
	struss(5) = concrete[5]->getStress();
	struss(6) = concrete[6]->getStress();
	sigl(0, 0) = struss(5);
	sigl(1, 1) = struss(6);
	tgtl(0, 0) = concrete[5]->getTangent();
	tgtl(1, 1) = -concrete[6]->getTangent();
	sigg = T * sigl * Ttran;
	tgtg = T * tgtl * Ttran;
	sig(1) += sigg(0, 0);
	sig(2) += sigg(1, 1);
	sig(4) += sigg(0, 1);
	K(1, 1) += tgtg(0, 0);
	K(2, 2) += tgtg(1, 1);
	K(4, 4) += tgtg(0, 1);

	epsl.Zero();
	sigl.Zero();
	tgtl.Zero();
	epsg(0, 0) = v(2) - (v(0) + v(2)) / 2;
	epsg(1, 1) = v(0) - (v(0) + v(2)) / 2;
	epsg(0, 1) = v(5);
	epsg(1, 0) = v(5);
	epsl = Ttran * epsg * T;
	etruss(7) = epsl(0, 0);
	etruss(8) = epsl(1, 1);
	concrete[7]->setTrialStrain(epsl(0, 0));
	concrete[8]->setTrialStrain(epsl(1, 1));
	struss(7) = concrete[7]->getStress();
	struss(8) = concrete[8]->getStress();
	sigl(0, 0) = struss(7);
	sigl(1, 1) = struss(8);
	tgtl(0, 0) = concrete[7]->getTangent();
	tgtl(1, 1) = -concrete[8]->getTangent();
	sigg = T * sigl * Ttran;
	tgtg = T * tgtl * Ttran;
	sig(2) += sigg(0, 0);
	sig(0) += sigg(1, 1);
	sig(5) += sigg(0, 1);
	K(2, 2) += tgtg(0, 0);
	K(0, 0) += tgtg(1, 1);
	K(5, 5) += tgtg(0, 1);

	return 0;
}

int
MCTM::setTrialStrain(const Vector& v, const Vector& rate)
{
	return setTrialStrain(v);
}

int
MCTM::setTrialStrainIncr(const Vector& v)
{
	opserr << "MCTM::setTrialStrainIncr -- subclass responsibility\n";
	exit(-1);
	return -1;
}

int
MCTM::setTrialStrainIncr(const Vector& v, const Vector& rate)
{
	opserr << "MCTM::setTrialStrainIncr -- subclass responsibility\n";
	exit(-1);
	return -1;
}

const Matrix&
MCTM::getTangent(void)
{
	return K;
}

const Matrix&
MCTM::getInitialTangent(void)
{
	if (Ki != 0)
		return *Ki;

	Ki = new Matrix(6, 6);

	(*Ki)(0, 0) = concrete[0]->getTangent();
	(*Ki)(1, 1) = concrete[1]->getTangent();
	(*Ki)(2, 2) = concrete[2]->getTangent();

	//Matrix tgtg(3, 3), tgtl(3, 3);
	//tgtl(0, 0) = concrete[3]->getInitialTangent();
	//tgtl(1, 1) = concrete[4]->getInitialTangent();
	//tgtg = T1 * tgtl * T1tran;
	//(*Ki) += tgtg;

	//tgtl.Zero();
	//tgtl(1, 1) = concrete[5]->getInitialTangent();
	//tgtl(2, 2) = concrete[6]->getInitialTangent();
	//tgtg = T2 * tgtl * T2tran;
	//(*Ki) += tgtg;

	//tgtl.Zero();
	//tgtl(0, 0) = concrete[7]->getInitialTangent();
	//tgtl(2, 2) = concrete[8]->getInitialTangent();
	//tgtg = T2 * tgtl * T2tran;
	//(*Ki) += tgtg;

	return (*Ki);
}

const Vector&
MCTM::getStress(void)
{
	return sig;
}

const Vector&
MCTM::getStrain(void)
{
	return eps;
}

int
MCTM::commitState(void)
{
	int res = 0;

	for (int i = 0; i < 9; i++)
		res += concrete[i]->commitState();

	return res;
}

int
MCTM::revertToLastCommit(void)
{
	int res = 0;

	for (int i = 0; i < 9; i++)
		res += concrete[i]->revertToLastCommit();

	return res;
}

int
MCTM::revertToStart(void)
{
	int res = 0;

	for (int i = 0; i < 9; i++)
		res += concrete[i]->revertToStart();

	return res;
}

NDMaterial*
MCTM::getCopy(void)
{
	opserr << "MCTM::getCopy -- subclass responsibility\n";
	exit(-1);
	return 0;
}

const char*
MCTM::getType(void) const
{
	return "MCTM";
}


Response*
MCTM::setResponse(const char** argv, int argc,
	OPS_Stream& output)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0) {
		const Vector& res = this->getStress();
		theResponse = new MaterialResponse(this, 1, this->getStress());
	}
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) {
		const Vector& res = this->getStrain();
		theResponse = new MaterialResponse(this, 2, this->getStress());
	}
	else if (strcmp(argv[0], "trussStress") == 0) {
		theResponse = new MaterialResponse(this, 3, struss);
	}
	else if (strcmp(argv[0], "trussStrain") == 0) {
		theResponse = new MaterialResponse(this, 4, etruss);
	}

	return theResponse;
}

int
MCTM::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case 1:
		return matInfo.setVector(this->getStress());
	case 2:
		return matInfo.setVector(this->getStrain());
	case 3:
		return matInfo.setVector(struss);
	case 4:
		return matInfo.setVector(etruss);


	default:
		return -1;
	}
}


int
MCTM::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int
MCTM::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
}

void
MCTM::Print(OPS_Stream& s, int flag)
{

}

int
MCTM::setParameter(const char** argv, int argc,
	Parameter& param)
{
	return 0;
}

int
MCTM::updateParameter(int parameterID, Information& info)
{
	return 0;
}

int
MCTM::activateParameter(int paramID)
{
	return 0;
}
