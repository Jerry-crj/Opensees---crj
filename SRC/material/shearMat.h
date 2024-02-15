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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/shearMat.h,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class definition for 
// shearMat. shearMat is based on an f2c of the FEDEAS material
// Concr2.f which is:
/*-----------------------------------------------------------------------
! concrete model with damage modulus
!       by MOHD YASSIN (1993)
! adapted to FEDEAS material library
! by D. Sze and Filip C. Filippou in 1994
-----------------------------------------------------------------------*/



#ifndef shearMat_h
#define shearMat_h

#include <Vector.h>
#include <Matrix.h>
#include <UniaxialMaterial.h>

#define oneOverThree 0.3
#define min(a,b) ( (a)<(b) ? (a):(b) )
#define max(a,b) ( (a)<(b) ? (b):(a) )

class shearMat : public UniaxialMaterial
{
public:
	shearMat(int tag, double fs0, double g, double gs, double perturb = 1.0e-7);

	shearMat(void);

	~shearMat();

	const char* getType(void) const { return "shearMat"; };

	UniaxialMaterial* getCopy(void);

	int setTrialStrain(double strain, double strainRate = 0.0);

	double getInitialTangent(void);

	double getStress(void);
	double getStrain(void);
	double getTangent(void);
//	double getRho(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int	getResponse(int responseID, Information& matInfo);

protected:

private:

	void getTrialStressTangent(double& st, double& tgt, double et);

	double fs0init;
	double epss0;
	double epssu;

	double fsp0P;
	double fsn0P;

	double ks;
	double kspP;
	double ksnP;
	double ksInit;

	double dmgsp;
	double dmgsn;
	double dmgspP;
	double dmgsnP;

	double eps;
	double sig;
	double epsP;
	double sigP;

	double epsressp;
	double epsressn;

	double Gs;

	double perturb;
};	


#endif
