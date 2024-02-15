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

// $Revision: 1.3 $
// $Date: 2007-06-08 00:38:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASCCTEST.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of ASCCTEST. 
// This ASCCTEST is based on an f2c of the FEDEAS material
// Concr2.f which is:
//-----------------------------------------------------------------------
// concrete model with damage modulus    
//       by MOHD YASSIN (1993)
// adapted to FEDEAS material library
// by D. Sze and Filip C. Filippou in 1994
//-----------------------------------------------------------------------


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ASCCTEST.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


#define MAT_TAG_ASCCTEST 123123141223

//Matrix ASCCTEST::e(2, 2);

void *
OPS_ASCCTEST()
{
	// Pointer to a uniaxial material that will be returned
	NDMaterial *theMaterial = 0;

	int iData[2];
	int numData = 2;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "Invalid #args, want: NDMaterial ASCCTEST " << iData[0] << " nDMaterial?\n";
		return 0;
	}
	
	NDMaterial *Mat2D = 0;

	Mat2D = OPS_getNDMaterial(iData[1]);

	// Parsing was successful, allocate the material
	theMaterial = new ASCCTEST(iData[0], Mat2D);

	if (theMaterial == 0) {
		opserr << "WARNING could not create NDMaterial of type ASCCTEST Material\n";
		return 0;
	}

	return theMaterial;
}

ASCCTEST::ASCCTEST(int tag, NDMaterial * _mat) :
	NDMaterial(tag, MAT_TAG_ASCCTEST), stress(3), strain(3), tangent(3, 3)
{
	theMat = _mat->getCopy();
	Matrix t(3, 3);
	t = _mat->getInitialTangent();
	tangent(0, 0) = 1.0e5;
	tangent(1, 1) = t(0, 0);
	tangent(1, 2) = t(0, 1);
	tangent(2, 1) = t(1, 0);
	tangent(2, 2) = t(1, 1);
}

ASCCTEST::ASCCTEST(void) :
	NDMaterial(0, MAT_TAG_ASCCTEST)
{

}

ASCCTEST::~ASCCTEST(void)
{
	// Does nothing
	delete theMat;
	theMat = 0;
}

NDMaterial*
ASCCTEST::getCopy(void)
{
	ASCCTEST *theCopy = new ASCCTEST(this->getTag(), theMat);

	return theCopy;
}

NDMaterial*
ASCCTEST::getCopy(const char *type)
{
	ASCCTEST *theCopy = new ASCCTEST(this->getTag(), theMat);

	return theCopy;
}

int
ASCCTEST::setTrialStrain(const Vector &trialStrain)
{
	strain = trialStrain;

	if (strain.Norm() > 0)
		bool debug = true;

	Vector tstrain(2), tstress(2);
	Matrix ttangent(2, 2);

	tstrain(0) = strain(1);
	tstrain(1) = strain(2);

	int err = theMat->setTrialStrain(tstrain);

	tstress = theMat->getStress();
	ttangent = theMat->getTangent();

	stress(0) = 0.0;
	stress(1) = tstress(0);
	stress(2) = tstress(1);

	tangent(0, 0) = 1.0e5;
	tangent(1, 1) = ttangent(0, 0);
	tangent(1, 2) = ttangent(0, 1);
	tangent(2, 1) = ttangent(1, 0);
	tangent(2, 2) = ttangent(1, 1);

	return err;
}

int 
ASCCTEST::setTrialStrain(const Vector &v, const Vector &r)
{
	return setTrialStrain(v);
}

int
ASCCTEST::setTrialStrainIncr(const Vector &v)
{
	Vector temp(3);
	temp = strain + v;
	return setTrialStrain(temp);
}

int
ASCCTEST::setTrialStrainIncr(const Vector &v, const Vector &rate)
{
	Vector temp(3);
	temp = strain + v;
	return setTrialStrain(temp);
}

const Vector&
ASCCTEST::getStrain(void)
{
	return strain;
}

const Vector&
ASCCTEST::getStress(void)
{
	return stress;
}

const Matrix&
ASCCTEST::getTangent(void)
{
	return tangent;
}

int
ASCCTEST::commitState(void)
{
	return theMat->commitState();
}

int
ASCCTEST::revertToLastCommit(void)
{
	return theMat->revertToLastCommit();
}

int
ASCCTEST::revertToStart(void)
{
	return theMat->revertToStart();
}

int
ASCCTEST::sendSelf(int commitTag, Channel &theChannel)
{
	//static Vector data(13);
	//data(0) = fc;
	//data(1) = epsc0;
	//data(2) = fcu;
	//data(3) = epscu;
	//data(4) = rat;
	//data(5) = ft;
	//data(6) = Ets;
	//data(7) = ecminP;
	//data(8) = deptP;
	//data(9) = epsP;
	//data(10) = sigP;
	//data(11) = eP;
	//data(12) = this->getTag();

	//if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
	//	opserr << "ASCCTEST::sendSelf() - failed to sendSelf\n";
	//	return -1;
	//}
	return 0;
}

int
ASCCTEST::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{

	//static Vector data(13);

	//if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
	//	opserr << "ASCCTEST::recvSelf() - failed to recvSelf\n";
	//	return -1;
	//}

	//fc = data(0);
	//epsc0 = data(1);
	//fcu = data(2);
	//epscu = data(3);
	//rat = data(4);
	//ft = data(5);
	//Ets = data(6);
	//ecminP = data(7);
	//deptP = data(8);
	//epsP = data(9);
	//sigP = data(10);
	//eP = data(11);
	//this->setTag(data(12));

	//e = eP;
	//sig = sigP;
	//eps = epsP;

	return 0;
}

void
ASCCTEST::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "ASCCTEST:(strain, stress, tangent) " << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ASCCTEST\", ";

	}
}



Response*
ASCCTEST::setResponse(const char **argv, int argc,
	OPS_Stream &output)
{
	Response *theResponse = 0;
	theResponse = theMat->setResponse(argv, argc, output);
	return theResponse;
}

int
ASCCTEST::getResponse(int responseID, Information &matInfo)
{
	return theMat->getResponse(responseID, matInfo);
}