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

// $Revision: 1.13 $
// $Date: 2006-09-05 21:21:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/MCTM.h,v $


#ifndef MCTM_h
#define MCTM_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for MCTMModel.
// MCTMModel is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) MCTM.h, revA"

#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#define MAT_TAG_MCTM 1945730

class MCTM : public NDMaterial
{
public:
	// Only called by subclasses to pass their tags to NDMaterialModel
	MCTM(int tag, UniaxialMaterial** c);
	MCTM(int tag, UniaxialMaterial** c, UniaxialMaterial** s);

	~MCTM(void);
	
	const char* getClassType(void) const { return "MCTM"; };
	
	double getRho();
	
	int setTrialStrain(const Vector& v);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);
	const Vector& getStress(void);
	const Vector& getStrain(void);
	
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// Create a copy of material parameters AND state variables
	// Called by GenericSectionXD
	NDMaterial* getCopy(void);

	// Create a copy of just the material parameters
	// Called by the continuum elements
	NDMaterial* getCopy(const char* type);

	// Return a string indicating the type of material model
	const char* getType(void) const;

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	Response* setResponse(const char** argv, int argc,
		OPS_Stream& s);
	int getResponse(int responseID, Information& matInformation);


	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	int activateParameter(int paramID);

protected:
	
	UniaxialMaterial** concrete;
	UniaxialMaterial** steel;

	double rho;

	Vector eps, sig, etruss, struss;
	Matrix K;
	Matrix* Ki;

	static Matrix T, Ttran;

private:

};


#endif
