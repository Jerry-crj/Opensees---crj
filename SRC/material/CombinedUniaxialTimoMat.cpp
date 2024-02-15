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

// $Revision: 1.5 $
// $Date: 2010-09-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/CombinedUniaxialTimoMat.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of CombinedUniaxialTimoMat. 
// This CombinedUniaxialTimoMat is based on an f2c of the FEDEAS material
// CombinedUniaxialTimoMat.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//-----------------------------------------------------------------------

#include <math.h>

#include <stdlib.h>
#include <CombinedUniaxialTimoMat.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <MaterialResponse.h>

# define MAT_TAG_CombinedUniaxialTimoMat 121253452

void *
OPS_CombinedUniaxialTimoMat()
{
	// Pointer to a uniaxial material that will be returned
	NDMaterial *theMaterial = 0;

	int iData[3];
	int numData = 3;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid nDMaterial CombinedUniaxialTimoMat tag" << endln;
		return 0;
	}

	UniaxialMaterial *mat1, *mat2;

	mat1 = OPS_GetUniaxialMaterial(iData[1]);
	mat2 = OPS_GetUniaxialMaterial(iData[2]);

	if (mat1 == 0 || mat2 == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type CombinedUniaxialTimoMat Material\n";
		opserr << "get UniaxialMaterials fialed\n";
		return 0;
	}

	UniaxialMaterial *theMat1, *theMat2;
	theMat1 = mat1->getCopy();
	theMat2 = mat2->getCopy();

		// Parsing was successful, allocate the material
	theMaterial = new CombinedUniaxialTimoMat(iData[0], theMat1, theMat2);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type CombinedUniaxialTimoMat Material\n";
		return 0;
	}

	return theMaterial;
}


CombinedUniaxialTimoMat::CombinedUniaxialTimoMat(int tag, UniaxialMaterial* mat1, UniaxialMaterial* mat2) :
	NDMaterial(tag, MAT_TAG_CombinedUniaxialTimoMat),
	sig(2), sigP(2), eps(2), epsP(2), k(2, 2), kP(2, 2), ki(0)
{
	theMaterials = new UniaxialMaterial * [2];
	theMaterials[0] = mat1;
	theMaterials[1] = mat2;

	k(0, 0) = mat1->getTangent();
	k(1, 1) = mat2->getTangent();

	kP = k;
}

CombinedUniaxialTimoMat::CombinedUniaxialTimoMat() :
	NDMaterial(0, MAT_TAG_CombinedUniaxialTimoMat),
	theMaterials(0), sig(2), sigP(2), eps(2), epsP(2), k(2, 2), kP(2, 2), ki(0)
{

}


CombinedUniaxialTimoMat::~CombinedUniaxialTimoMat(void)
{
	for (int i = 0; i < 2; i++)
		delete theMaterials[i];

	delete theMaterials;
}

NDMaterial*
CombinedUniaxialTimoMat::getCopy(void)
{
	UniaxialMaterial* newMat1, * newMat2;
	newMat1 = theMaterials[0]->getCopy();
	newMat2 = theMaterials[1]->getCopy();
	CombinedUniaxialTimoMat *theCopy = new CombinedUniaxialTimoMat(this->getTag(), newMat1, newMat2);

	return theCopy;
}

NDMaterial*
CombinedUniaxialTimoMat::getCopy(const char *type)
{
	NDMaterial* theCopy = this->getCopy();
	return theCopy;
}

const Matrix&
CombinedUniaxialTimoMat::getInitialTangent(void)
{
	if (ki != 0)
		return *ki;

	(*ki)(0, 0) = theMaterials[0]->getInitialTangent();
	(*ki)(1, 1) = theMaterials[0]->getInitialTangent();
	
	return *ki;
}

int
CombinedUniaxialTimoMat::setTrialStrain(const Vector &trialStrian)
{
	eps = trialStrian;

	int err = 0;
	err += theMaterials[0]->setTrialStrain(eps(0));
	err += theMaterials[1]->setTrialStrain(eps(1));

	sig(0) = theMaterials[0]->getStress();
	sig(1) = theMaterials[1]->getStress();

	k(0, 0) = theMaterials[0]->getTangent();
	k(1, 1) = theMaterials[1]->getTangent();

	return err;
}



const Vector&
CombinedUniaxialTimoMat::getStrain(void)
{
	return eps;
}

const Vector&
CombinedUniaxialTimoMat::getStress(void)
{
	return sig;
}

const Matrix&
CombinedUniaxialTimoMat::getTangent(void)
{
	return k;
}

int
CombinedUniaxialTimoMat::commitState(void)
{
	for (int i = 0; i < 2; i++)
		theMaterials[i]->commitState();

	epsP = eps;
	sigP = sig;
	kP = k;
	
	return 0;
}

int
CombinedUniaxialTimoMat::revertToLastCommit(void)
{
	for (int i = 0; i < 2; i++)
		theMaterials[i]->revertToLastCommit();

	eps = epsP;
	sig = sigP;
	k = kP;

	return 0;
}

int
CombinedUniaxialTimoMat::revertToStart(void)
{
	k.Zero();
	k(0, 0) = theMaterials[0]->getTangent();
	k(1, 1) = theMaterials[1]->getTangent();

	kP = k;

	sig.Zero();
	eps.Zero();

	sigP.Zero();
	epsP.Zero();

	return 0;
}

int
CombinedUniaxialTimoMat::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}


int
CombinedUniaxialTimoMat::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
}


int
CombinedUniaxialTimoMat::getResponse(int responseID, Information& matInfo)
{
	switch (responseID)
	{
	case 1:
	{
		Vector res(4);
		res(0) = sig(0);
		res(1) = sig(1);
		res(2) = eps(0);
		res(3) = eps(1);
		return matInfo.setVector(res);
	}
	default:
		return -1;
	}
}


Response*
CombinedUniaxialTimoMat::setResponse(const char** argv, int argc,
	OPS_Stream& output)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "info") == 0) {
		Vector res(4);
		theResponse = new MaterialResponse(this, 1, res);
	}

	return theResponse;
}


void
CombinedUniaxialTimoMat::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		//    s << "CombinedUniaxialTimoMat:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
		s << "\tMaterial type: CombinedUniaxialTimoMat, Material tag: " << this->getTag() << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"CombinedUniaxialTimoMat\", ";
	}
}


UniaxialMaterial**
CombinedUniaxialTimoMat::getMaterials(void)
{
	return theMaterials;
}