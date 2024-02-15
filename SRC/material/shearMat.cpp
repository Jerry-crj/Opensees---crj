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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/shearMat.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of shearMat. 
// This shearMat is based on an f2c of the FEDEAS material
// Concr2.f which is:
//-----------------------------------------------------------------------
// concrete model with damage modulus    
//       by MOHD YASSIN (1993)
// adapted to FEDEAS material library
// by D. Sze and Filip C. Filippou in 1994
//-----------------------------------------------------------------------


#include <stdlib.h>
#include <string>
#include <math.h>

#include <shearMat.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <MaterialResponse.h>

#define MAT_TAG_shearMat 12312314123211

void*
OPS_shearMat()
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	int    iData[1];
	double dData[4];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid UniaxialMaterial shearMat tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData < 3) {
		opserr << "Invalid #args, want: UniaxialMaterial shearMat " << iData[0] << " fpc? epsc0? epscu? epst0? epstu?\n";
		return 0;
	}

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid #args, want: UniaxialMaterial shearMat " << iData[0] << " fpc? epsc0? epscu? epst0? epstu?\n";
		return 0;
	}

	// Parsing was successful, allocate the material
	if (numData == 3)
	{
		theMaterial = new shearMat(iData[0], dData[0], dData[1], dData[2]);
	}

	if (theMaterial == 0) {
		opserr << "WARNING could not create UniaxialMaterial of type shearMat Material\n";
		return 0;
	}

	return theMaterial;
}

shearMat::shearMat(int tag, double ss0, double g, double gs, double perturb) :
	UniaxialMaterial(tag, MAT_TAG_shearMat), fs0init(ss0),
	dmgsp(0), dmgspP(0), dmgsn(0), dmgsnP(0),
	epsressp(0), epsressn(0), eps(0), epsP(0), sig(0), sigP(0), perturb(1e-7)
{
	Gs = gs;

	ksInit = g;
	kspP = g;
	ksnP = g;
	ks = g;

	epss0 = abs(ss0 / ksInit);
	epssu = epss0 + abs(ss0 / gs);

	fsp0P = abs(ss0);
	fsn0P = -abs(ss0);

	epsressp = epss0;
	epsressn = -epss0;
}

shearMat::shearMat(void) :
	UniaxialMaterial(0, MAT_TAG_shearMat),
	epss0(0), epssu(0), fs0init(0),
	fsp0P(0), fsn0P(0), dmgsp(0), dmgspP(0), dmgsn(0), dmgsnP(0),
	ksInit(0), epsressp(0), epsressn(0),
	Gs(0), epsP(2), sigP(2), eps(2), sig(2), ks(0), ksnP(0), kspP(0), 
	perturb(1e-7)
{

}

shearMat::~shearMat(void)
{
	// Does nothing
}

UniaxialMaterial*
shearMat::getCopy(void)
{
	shearMat* theCopy = new shearMat(this->getTag(), fs0init, ksInit, Gs, perturb);

	return theCopy;
}

int
shearMat::setTrialStrain(double strain, double strainRate)
{
	eps = strain;
	sig = sigP + ksInit * (eps - epsP);
	ks = ksInit;

	if (eps > 0)
	{
		if (eps > epsressp)
		{
			if (eps > epssu)
			{
				sig = 0;
				ks = 1.0e-10;
			}
			else
			{
				sig = fs0init * (1 - (eps - epss0) / (epssu - epss0));
				ks = -Gs;
			}
		}
		else
		{
			double sigMax = kspP * eps;
			if (sig > sigMax)
			{
				sig = sigMax;
				ks = kspP;
			}
			if (sig * eps < 0)
			{
				sig = 0.0;
				ks = 1.0e-10;
			}
		}
	}
	else
	{
		if (eps < epsressn)
		{
			if (eps < -epssu)
			{
				sig = 0;
				ks = 1.0e-10;
			}
			else
			{
				sig = -fs0init * (1 - abs(eps + epss0) / (epssu - epss0));
				ks = -Gs;
			}
		}
		else
		{
			double sigMax = ksnP * eps;
			if (sig < sigMax)
			{
				sig = sigMax;
				ks = ksnP;
			}
			if (sig * eps < 0)
			{
				sig = 0.0;
				ks = 1.0e-10;
			}
		}
	}

	return 0;
}

double
shearMat::getStrain(void)
{
	return eps;
}

double
shearMat::getStress(void)
{
	return sig;
}

double
shearMat::getTangent(void)
{
	return ks;
}

double
shearMat::getInitialTangent(void)
{
	return ksInit;
}


int
shearMat::commitState(void)
{
	if (eps > epsressp)
	{
		epsressp = eps;
		kspP = sig / eps;
	}
	else if (eps < epsressn)
	{
		epsressn = eps;
		ksnP = sig / eps;
	}

	dmgsnP = dmgsn;
	dmgspP = dmgsp;

	sigP = sig;
	epsP = eps;


	return 0;
}

int
shearMat::revertToLastCommit(void)
{
	dmgsn = dmgsnP;
	dmgsp = dmgspP;

	sig = sigP;
	eps = epsP;

	return 0;
}

int
shearMat::revertToStart(void)
{

	dmgsp = 0.0;
	dmgsn = 0.0;
	dmgspP = 0.0;
	dmgsnP = 0.0;

	fsp0P = fs0init;
	fsn0P = -fs0init;

	sig = 0.0;
	eps = 0.0;
	sigP = 0.0;
	epsP = 0.0;

	return 0;
}

int
shearMat::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int
shearMat::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
}

void
shearMat::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "shearMat:(strain, stress, tangent) " << eps << " " << sig << " " << ks << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"shearMat\", ";
	}
}

Response*
shearMat::setResponse(const char** argv, int argc,
	OPS_Stream& output)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0) {
		return theResponse = this->UniaxialMaterial::setResponse(argv, argc, output);
	}
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) {
		return theResponse = this->UniaxialMaterial::setResponse(argv, argc, output);
	}
	if (strcmp(argv[0], "dmgfct") == 0) {
		Vector res(6);
		theResponse = new MaterialResponse(this, 3, res);
	}
	else if (strcmp(argv[0], "tangent") == 0) {
		Vector res(2);
		theResponse = new MaterialResponse(this, 4, res);
	}

	return theResponse;
}

int
shearMat::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case 1:
		return matInfo.setVector(this->getStress());
	case 2:
		return matInfo.setVector(this->getStrain());
	case 3:
	{
		Vector res(6);
		return matInfo.setVector(res);
	}
	case 4:
	{
		Vector res(2);
		return matInfo.setVector(res);
	}
	default:
		return -1;
	}
}


void
shearMat::getTrialStressTangent(double& st, double& tgt, double et)
{
	st = sigP + ksInit * (et - epsP);
	ks = ksInit;

	if (st > 0)
	{
		if (st > epsressp)
		{
			sig = fs0init * (et - epss0) / (epssu - epss0);
			ks = -Gs;
		}
		else
		{
			double sigMax = kspP * et;
			if (st > sigMax)
			{
				st = sigMax;
				ks = kspP;
			}
			if (st * et < 0)
			{
				st = 0.0;
				ks = 1.0e-10;
			}
		}
	}
	else
	{
		if (st < epsressn)
		{
			sig = -fs0init * abs(et + epss0) / (epssu - epss0);
			ks = -Gs;
		}
		else
		{
			double sigMax = ksnP * et;
			if (st < sigMax)
			{
				st = sigMax;
				ks = ksnP;
			}
			if (st * et < 0)
			{
				st = 0.0;
				ks = 1.0e-10;
			}
		}
	}
}