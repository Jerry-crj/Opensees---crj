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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ovalConcrete.h,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class definition for 
// ovalConcrete. ovalConcrete is based on an f2c of the FEDEAS material
// Concr2.f which is:
/*-----------------------------------------------------------------------
! concrete model with damage modulus
!       by MOHD YASSIN (1993)
! adapted to FEDEAS material library
! by D. Sze and Filip C. Filippou in 1994
-----------------------------------------------------------------------*/



#ifndef ovalConcrete_h
#define ovalConcrete_h

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>

#define oneOverThree 0.3
#define min(a,b) ( (a)<(b) ? (a):(b) )
#define max(a,b) ( (a)<(b) ? (b):(a) )

class ovalConcrete : public NDMaterial
{
public:
	ovalConcrete(int tag, double fc, double ec0, double fcu, double ecu,
		double rat, double ft, double ets, double fs0, double grat, double gs,
		double rmax = 0.3, double rmin = 0.1, double p = 1.0e-7);

	ovalConcrete(void);

	~ovalConcrete();

	const char* getType(void) const { return "ovalConcrete"; };
	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* type);

	int setTrialStrain(const Vector& v);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	const Vector& getStress(void);
	const Vector& getStrain(void);
	double getRho(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);


	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int	getResponse(int responseID, Information& matInfo);
	//void setPerturbValue(double pv) { perturb = pv; }
	//Vector getDmgage(void) { Vector res(2); res(0) = dmgtP; res(1) = dmgcP; return res; }

protected:

private:

	void getDltDmgTensile(const Vector line0, const Vector lineu, Vector sig, double& dmgt, double& dmgs);
	double getDltDmgShear(const Vector line0, const Vector lineu, Vector sig);
	double getDltDmgCompressive(const Vector ellipse0, const Vector ellipseu, Vector sig);

	void getStressTensile(Vector sig, Vector& sigCor, double dmgt, double dmgs);
	void getStressShear(Vector sig, Vector& sigCor, double dmgc, double dmgt, double dmgs);
	void getStressCompressive(Vector sig, Vector& sigCor, double dmgc, double dmgs);

	Matrix getTangentTensile(Vector sig, Vector sigCor, const Vector line0P, const Vector lineuP);
	Matrix getTangentShear(Vector sig, Vector sigCor, const Vector line0P, const Vector lineuP);
	Matrix getTangentShear(Vector sig, Vector sigCor, const Vector line0P);
	Matrix getTangentCopressive(Vector sig, Vector sigCor, const Vector ellipse0P, const Vector ellipseuP);
	Matrix getTangentCopressive(Vector sig, Vector sigCor, const Vector ellipse0P);

	//////////////////   deal geometry method
	const Vector getVerticalLine(Vector v1);
	const Vector getLineByTwoPoints(Vector v1, Vector v2);
	void getEllipseLineCrossPoint(Vector& crossP1, Vector& crossP2, Vector ellipse, Vector line);
	void getTwoLinesCrossPoint(Vector& crossP, Vector line1, Vector line2);

	void getTrialStressTangent(Vector& str, Matrix& tgt, Vector etr);

	double fcuinit;
	double fc0init;
	double ft0init;
	double fs0init;
	double epsc0;
	double epscu;
	double epst0;
	double epstu;
	double epss0;
	double epssu;

	double ft0;
	double fs0;
	double fc0;

	double ft0P;
	double fsp0P;
	double fsn0P;
	double fc0P;
	double ftiP;
	double fciP;
	double ftspiP;
	double ftsniP;
	double fcspiP;
	double fcsniP;

	double ktP;
	double kcP;
	double ktspP;
	double ktsnP;
	double kcspP;
	double kcsnP;
	double ksInit;
	double knInit;

	double dmgc;
	double dmgt;
	double dmgcsp;
	double dmgcsn;
	double dmgtsp;
	double dmgtsn;

	double dmgcP;
	double dmgtP;
	double dmgcspP;
	double dmgcsnP;
	double dmgtspP;
	double dmgtsnP;

	Vector eps;
	Vector sig;
	Vector epsP;
	Vector sigP;

	Vector sigCor;

	Matrix e;
	Matrix e0;
	Matrix eP;

	double perturb;

	double unldrat;
	double ultmDmg;

	double axRatpP;
	double axRatnP;
	double axRatMax;
	double axRatMin;

	double epsp;
	double sigp;
	double epsresn;
	double epsressp;
	double epsressn;
	double epsresnMin;
	double epsresspMin;
	double epsressnMin;

	double Ets;
	double Gs;

	double rho;
};	


#endif
