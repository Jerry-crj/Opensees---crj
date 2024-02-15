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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/QuadToTimoMat.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of QuadToTimoMat. 
// This QuadToTimoMat is based on an f2c of the FEDEAS material
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

#include <QuadToTimoMat.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


#define MAT_TAG_QuadToTimoMat 12312314122312

//Matrix QuadToTimoMat::e(2, 2);

void *
OPS_QuadToTimoMat()
{
	// Pointer to a uniaxial material that will be returned
	NDMaterial *theMaterial = 0;

	int iData[2];
	int numData = 2;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "Invalid #args, want: NDMaterial QuadToTimoMat " << iData[0] << " nDMaterial?\n";
		return 0;
	}
	
	NDMaterial *Mat2D = 0;

	Mat2D = OPS_getNDMaterial(iData[1]);

	// Parsing was successful, allocate the material
	theMaterial = new QuadToTimoMat(iData[0], Mat2D);

	if (theMaterial == 0) {
		opserr << "WARNING could not create NDMaterial of type QuadToTimoMat Material\n";
		return 0;
	}

	return theMaterial;
}

QuadToTimoMat::QuadToTimoMat(int tag, NDMaterial * quadMat) :
	NDMaterial(tag, MAT_TAG_QuadToTimoMat), stress(2), strain(2), tangent(2, 2)
{
	//theMat = (ElasticIsotropicPlaneStress2D *) quadMat->getCopy("PlaneStress2D");
	theMat = (ElasticIsotropicThreeDimensional *) quadMat->getCopy("3D");
	Matrix t = theMat->getInitialTangent();
	tangent(0, 0) = t(2, 2);
	tangent(0, 1) = t(2, 3);
	tangent(1, 0) = t(3, 2);
	tangent(1, 1) = t(3, 3);
}

QuadToTimoMat::QuadToTimoMat(void) :
	NDMaterial(0, MAT_TAG_QuadToTimoMat)
{

}

QuadToTimoMat::~QuadToTimoMat(void)
{
	// Does nothing
	delete theMat;
	theMat = 0;
}

NDMaterial*
QuadToTimoMat::getCopy(void)
{
	QuadToTimoMat *theCopy = new QuadToTimoMat(this->getTag(), theMat);

	return theCopy;
}

NDMaterial*
QuadToTimoMat::getCopy(const char *type)
{
	QuadToTimoMat *theCopy = new QuadToTimoMat(this->getTag(), theMat);

	return theCopy;
}

int
QuadToTimoMat::setTrialStrain(const Vector &trialStrain)
{
	strain = trialStrain;

	if (strain.Norm() > 0)
		bool debug = true;

	Vector tstrain(6), tstress(6);
	Matrix ttangent(6, 6);

	tstrain(2) = strain(0);
	tstrain(3) = strain(1);

	int err = theMat->setTrialStrain(tstrain);

	tstress = theMat->getStress();
	ttangent = theMat->getTangent();

	stress(0) = tstress(2);
	stress(1) = tstress(3);

	tangent(0, 0) = ttangent(2, 2);
	tangent(0, 1) = ttangent(2, 3);
	tangent(1, 0) = ttangent(3, 2);
	tangent(1, 1) = ttangent(3, 3);

	return err;
}

int 
QuadToTimoMat::setTrialStrain(const Vector &v, const Vector &r)
{
	return setTrialStrain(v);
}

int
QuadToTimoMat::setTrialStrainIncr(const Vector &v)
{
	Vector temp(3);
	temp = strain + v;
	return setTrialStrain(temp);
}

int
QuadToTimoMat::setTrialStrainIncr(const Vector &v, const Vector &rate)
{
	Vector temp(3);
	temp = strain + v;
	return setTrialStrain(temp);
}

const Vector&
QuadToTimoMat::getStrain(void)
{
	return strain;
}

const Vector&
QuadToTimoMat::getStress(void)
{
	return stress;
}

const Matrix&
QuadToTimoMat::getTangent(void)
{
	return tangent;
}

int
QuadToTimoMat::commitState(void)
{
	return theMat->commitState();
}

int
QuadToTimoMat::revertToLastCommit(void)
{
	return theMat->revertToLastCommit();
}

int
QuadToTimoMat::revertToStart(void)
{
	return theMat->revertToStart();
}

int
QuadToTimoMat::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int
QuadToTimoMat::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	return 0;
}

void
QuadToTimoMat::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "\t\tQuadToTimoMat: " << endln;
		s << "\t\tstrain:  " << strain(0) << "\t" << strain(1) << endln;
		s << "\t\tstress:  " << stress(0) << "\t" << stress(1) << endln;
		s << "\t\ttangent:  " << tangent(0, 0) << "\t" << tangent(0, 1) << endln;
		s << "\t\ttangent:  " << tangent(1, 0) << "\t" << tangent(1, 1) << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"TimoToQuadMat\", ";
	}
}

Response*
QuadToTimoMat::setResponse(const char **argv, int argc,
	OPS_Stream &output)
{
	Response *theResponse = 0;
	theResponse = theMat->setResponse(argv, argc, output);
	return theResponse;
}

int
QuadToTimoMat::getResponse(int responseID, Information &matInfo)
{
	return theMat->getResponse(responseID, matInfo);
}