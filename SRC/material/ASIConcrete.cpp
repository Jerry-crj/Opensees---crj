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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASIConcrete.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of ASIConcrete. 
// This ASIConcrete is based on an f2c of the FEDEAS material
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

#include <ASIConcrete.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <MaterialResponse.h>

#define MAT_TAG_ASIConcrete 123123141232

void*
OPS_ASIConcrete()
{
	// Pointer to a uniaxial material that will be returned
	NDMaterial* theMaterial = 0;

	int    iData[1];
	double dData[13];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid NDMaterial ASIConcrete tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData < 10) {
		opserr << "Invalid #args, want: NDMaterial ASIConcrete " << iData[0] << " fpc? epsc0? epscu? epst0? epstu?\n";
		return 0;
	}

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid #args, want: NDMaterial ASIConcrete " << iData[0] << " fpc? epsc0? epscu? epst0? epstu?\n";
		return 0;
	}

	// Parsing was successful, allocate the material
	if (numData == 10)
	{
		theMaterial = new ASIConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9]);
	}
	else if (numData == 11)
	{
		theMaterial = new ASIConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9], dData[10]);
	}
	else if (numData == 12)
	{
		theMaterial = new ASIConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);
	}
	else if (numData == 13)
	{
		theMaterial = new ASIConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]);
	}

	if (theMaterial == 0) {
		opserr << "WARNING could not create NDMaterial of type ASIConcrete Material\n";
		return 0;
	}

	return theMaterial;
}

ASIConcrete::ASIConcrete(int tag, double sc0, double ec0, double scu, double ecu, double rat,
	double st, double ets, double ss0, double grat, double gs, double rmax, double rmin, double p) :
	NDMaterial(tag, MAT_TAG_ASIConcrete), dmgt(0), dmgtP(0), dmgc(0), dmgcP(0), dmgs(0), dmgsP(0),
	epsresn(0), epsresnMin(0), epsressp(0), epsresspMin(0), epsressn(0), epsressnMin(0),
	eps(2), epsP(2), sig(2), sigCor(2), sigP(2), e(2, 2), e0(2, 2), eP(2, 2), lambda(rat), rho(0)
{
	Ets = ets;
	Gs = gs;

	ft = abs(st);
	fs = abs(ss0);
	fc = -abs(sc0);
	fcu = -abs(scu);

	kcP = fc / ec0;
	ktP = kcP;
	ktsP = kcP * grat;
	kcsP = kcP * grat;
	knInit = kcP;
	ksInit = ktsP;

	epsc0 = -abs(ec0);
	epscu = -abs(ecu);
	epst0 = abs(st / kcP);
	epstu = epst0 + abs(st / Ets);
	epss0 = abs(ss0 / ktsP);
	epssu = epss0 + abs(ss0 / gs);

	fty = ft;
	fcy = fc;
	fsy = fs;

	ftyP = abs(fty);
	fsyP = abs(fsy);
	fcyP = -abs(fcy);
	ftiP = ftyP + ktP * (epstu - epst0);
	fciP = fcyP + kcP * (epscu - epsc0);
	fsiP = fsyP + ktsP * (epssu - epss0);

	ultmDmg = 1 - fcu / fc;

	eP(0, 0) = ktP;
	eP(1, 1) = ktsP;

	epsp = (ecu * kcP * lambda - scu) / (kcP * (-1 + lambda));
	sigp = (ecu * kcP * lambda - scu) / (-1 + lambda);

	perturb = p;

	aspectRatio  = rmax;
	aspectRatioInit = rmax;
	aspectRatioUlti = rmin;
}

ASIConcrete::ASIConcrete(void) :
	NDMaterial(0, MAT_TAG_ASIConcrete),
	fcy(0), epsc0(0), epscu(0), epst0(0), epstu(0), epss0(0), epssu(0),
	fc(0), ft(0), fs(0), fcu(0), fty(0), fsy(0),
	fcyP(0), ftyP(0), fsyP(0), fciP(0), ftiP(0), fsiP(0), 
	dmgt(0), dmgtP(0), dmgc(0), dmgcP(0), dmgs(0), dmgsP(0),
	ktsP(0), kcsP(0), kcP(0), ktP(0), ksInit(0), knInit(0),
	epsresn(0), epsresnMin(0), epsressp(0), epsresspMin(0), epsressn(0), epsressnMin(0),
	Ets(0), Gs(0), epsP(2), sigP(2), eps(2), sig(2), e(2, 2), eP(2, 2),
	perturb(1e-7), aspectRatio (0.3), aspectRatioInit(0.3), aspectRatioUlti(0.3),
	epsp(0.0), sigp(0.0), lambda(0), ultmDmg(1.0), rho(0)
{

}

ASIConcrete::~ASIConcrete(void)
{
	// Does nothing
}

NDMaterial*
ASIConcrete::getCopy(void)
{
	ASIConcrete* theCopy = new ASIConcrete(this->getTag(), fc, epsc0, fcu, epscu,
		lambda, ft, Ets, fs, ktsP / ktP, Gs, aspectRatioInit, aspectRatioUlti, perturb);

	return theCopy;
}

NDMaterial*
ASIConcrete::getCopy(const char* type)
{
	ASIConcrete* theCopy = 0;
	theCopy = (ASIConcrete*)this->getCopy();

	return theCopy;
}


int
ASIConcrete::setTrialStrain(const Vector& trialStrain)
{
	eps = trialStrain;

	getTrialStressTangent(sig, e, eps);

	if (dmgtP > 1.0 - 1.0e-10 && sig(0) > 0)
	{
		dmgt = 1.0;
		sig(0) = 0;
	}

	Vector sigCor(2);

	if (eps(0) >= epsresn)
	{
		if (dmgt == 1.0)
		{
			sig.Zero();
			e.Zero();
			e(0, 0) = 1.0e-10;
			e(1, 1) = 1.0e-10;
		}
		else
		{
			if (sig(0) / ftyP + abs(sig(1) / fsyP) - 1 > 0.0)
			{
				Vector p1(2), p2(2), line0P(3), lineiP(3);

				if (sig(1) > 0)
				{
					p1(0) = ftyP;
					p2(1) = fsyP;
					line0P = getLineByTwoPoints(p1, p2);
					p1(0) = ftiP;
					p2(1) = fsiP;
					lineiP = getLineByTwoPoints(p1, p2);
				}
				else
				{
					p1(0) = ftyP;
					p2(1) = -fsyP;
					line0P = getLineByTwoPoints(p1, p2);
					p1(0) = ftiP;
					p2(1) = -fsiP;
					lineiP = getLineByTwoPoints(p1, p2);
				}

				getDmgTensile(line0P, lineiP, sig, dmgt);

				getStressTensile(sig, sigCor, dmgt);

				e = getTangentTensile(sig, sigCor, line0P, lineiP);

				sig = sigCor;
			}
		}
	}
	else if (sig(0) > fcyP / 2 && eps(0) < epsresn)
	{
		Vector p1(2), p2(2), line0P(3), lineiP(3);

		if (sig(1) > 0)
		{
			p1(0) = fcyP / 2;
			p1(1) = abs(fcyP) / 2 * aspectRatio ;
			p2(1) = fsyP;
			line0P = getLineByTwoPoints(p1, p2);
		}
		else
		{
			p1(0) = fcyP / 2;
			p1(1) = -abs(fcyP) / 2 * aspectRatio ;
			p2(1) = -fsyP;
			line0P = getLineByTwoPoints(p1, p2);
		}
		
		if ((line0P(0) * sig(0) + line0P(1) * sig(1) + line0P(2) > 0 && sig(1) > 0)
			|| (line0P(0) * sig(0) + line0P(1) * sig(1) + line0P(2) < 0 && sig(1) < 0))// out of yield surface
		{
			if (dmgtP > 1.0 - 1.0e-10 && dmgsP > ultmDmg - 1.0e-10)
			{
				dmgs = 1.0;
				dmgs = ultmDmg;

				getStressShear(sig, sigCor, dmgcP, 1.0, ultmDmg);

				e = getTangentShear(sig, sigCor, line0P);

				sig = sigCor;
			}
			else
			{
				if (sig(1) > 0)
				{
					p1(0) = fcyP / 2;
					p1(1) = abs(fciP - fcyP / 2) * aspectRatio ;
					p2(1) = fsiP;
					lineiP = getLineByTwoPoints(p1, p2);
				}
				else
				{
					p1(0) = fcyP / 2;
					p1(1) = -abs(fciP - fcyP / 2) * aspectRatio ;
					p2(1) = -fsiP;
					lineiP = getLineByTwoPoints(p1, p2);
				}

				getDmgShear(line0P, lineiP, sig, dmgt, dmgs);

				getStressShear(sig, sigCor, dmgcP, dmgt, dmgs);

				e = getTangentShear(sig, sigCor, line0P, lineiP);

				sig = sigCor;
			}
		}
	}
	else   // shear and compressive area
	{
		Vector ellipse0P(3);
		ellipse0P(0) = abs(fcyP) / 2;
		ellipse0P(1) = ellipse0P(0) * aspectRatio ;
		ellipse0P(2) = fcyP / 2;

		if (pow((sig(0) - ellipse0P(2)) / ellipse0P(0), 2) + pow(sig(1) / ellipse0P(1), 2) - 1 > 0.0)  // outside of the ultimate sufrace
		{
			if (dmgcP > ultmDmg - 1.0e-10 && dmgsP > ultmDmg - 1.0e-10)
			{
				getStressCompressive(sig, sigCor, ultmDmg, ultmDmg);

				e = getTangentCopressive(sig, sigCor, ellipse0P);

				sig = sigCor;
			}
			else
			{
				Vector ellipseiP(3);
				ellipseiP(0) = abs(fciP - fcyP / 2);
				ellipseiP(1) = ellipseiP(0) * aspectRatio ;
				ellipseiP(2) = fcyP / 2;

				getDmgCompressive(ellipse0P, ellipseiP, sig, dmgc, dmgs);

				getStressCompressive(sig, sigCor, dmgc, dmgs);

				e = getTangentCopressive(sig, sigCor, ellipse0P, ellipseiP);

				sig = sigCor;
			}
		}

	}

	return 0;
}


int
ASIConcrete::setTrialStrain(const Vector& v, const Vector& r)
{
	return setTrialStrain(v);
}

int
ASIConcrete::setTrialStrainIncr(const Vector& v)
{
	Vector temp(3);
	temp = eps + v;
	return setTrialStrain(temp);
}

int
ASIConcrete::setTrialStrainIncr(const Vector& v, const Vector& rate)
{
	Vector temp(3);
	temp = eps + v;
	return setTrialStrain(temp);
}

const Vector&
ASIConcrete::getStrain(void)
{
	return eps;
}

const Vector&
ASIConcrete::getStress(void)
{
	return sig;
}

const Matrix&
ASIConcrete::getTangent(void)
{
	return e;
}

double
ASIConcrete::getRho(void)
{
	return rho;
}

const Matrix&
ASIConcrete::getInitialTangent(void)
{
	return e0;
}


int
ASIConcrete::commitState(void)
{
	sigP = sig;
	epsP = eps;

	aspectRatio  = (aspectRatioUlti - aspectRatioInit) / ultmDmg * dmgs + aspectRatioInit;

	ftyP = ft * (1 - dmgt);
	fcyP = fc * (1 - dmgc);
	fsyP = fs * (1 - dmgt);

	ktP = ft * (1 - dmgt) / (epst0 + dmgt * (epstu - epst0));
	kcP = (fc * (1 - dmgc) - sigp) / (epsc0 + dmgc * (epscu - epsc0) / ultmDmg - epsp);
	kcsP = fs * (1 - dmgs) / (epss0 + dmgs * (epssu - epss0));
	ktsP = fs * (1 - dmgt) / (epss0 + dmgt * (epssu - epss0));


	ftiP = ftyP + ktP * (epstu - epst0) * (1 - dmgt);
	fciP = fcyP + kcP * (epscu - epsc0) / ultmDmg * (1 - dmgc);
	fsiP = fsyP + kcsP * (epssu - epss0) * (1 - dmgt);

	if (dmgc > dmgcP)
	{
		if (epsresnMin > eps(0) - sig(0) / kcP)
			epsresnMin = eps(0) - sig(0) / kcP;
	}
	else if (dmgc == ultmDmg)
	{
		if (sig(0) < fcyP / 2 &&
			pow((sig(0) - fcyP / 2) / abs(fcyP / 2), 2) + pow(sig(1) / abs(fcyP / 2) / aspectRatio , 2) - 1 > -1.0e-10)
			if (epsresnMin > eps(0) - sig(0) / kcP)
				epsresnMin = eps(0) - sig(0) / kcP;
	}
	epsresn = epsresnMin;

	if (eps(0) > epsresn)
	{
		eP(0, 0) = ktP;
		eP(1, 1) = ktsP;
	}
	else
	{
		eP(0, 0) = kcP;
		eP(1, 1) = kcsP;
	}

	dmgtP = dmgt;
	dmgcP = dmgc;
	dmgsP = dmgs;

	return 0;
}

int
ASIConcrete::revertToLastCommit(void)
{
	dmgt = dmgtP;
	dmgc = dmgcP;
	dmgs = dmgsP;

	sig = sigP;
	eps = epsP;

	return 0;
}

int
ASIConcrete::revertToStart(void)
{
	dmgt = 0.0;
	dmgc = 0.0;
	dmgs = 0.0;
	dmgtP = 0.0;
	dmgcP = 0.0;
	dmgsP = 0.0;

	ftyP = ft;
	fsyP = fs;
	fcyP = fc;
	ftiP = ftyP + ktP * (epstu - epst0);
	fciP = fcyP + kcP * (epscu - epsc0);
	fsiP = fsyP + kcsP * (epssu - epss0);

	sig.Zero();
	eps.Zero();
	sigP.Zero();
	epsP.Zero();

	return 0;
}

int
ASIConcrete::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int
ASIConcrete::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
}

void
ASIConcrete::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "ASIConcrete:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ASIConcrete\", ";
		s << "\"Ec\": " << fcy / epsc0 << ", ";
		s << "\"fc\": " << fcy << ", ";
		s << "\"epsc\": " << epsc0 << ", ";
		s << "\"epscu\": " << epscu << ", ";
	}
}

Response*
ASIConcrete::setResponse(const char** argv, int argc,
	OPS_Stream& output)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0) {
		return theResponse = this->NDMaterial::setResponse(argv, argc, output);
	}
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0) {
		return theResponse = this->NDMaterial::setResponse(argv, argc, output);
	}
	if (strcmp(argv[0], "dmgfct") == 0) {
		Vector res(3);
		theResponse = new MaterialResponse(this, 3, res);
	}
	else if (strcmp(argv[0], "tangent") == 0) {
		Vector res(4);
		theResponse = new MaterialResponse(this, 4, res);
	}

	return theResponse;
}

int
ASIConcrete::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case 1:
		return matInfo.setVector(this->getStress());
	case 2:
		return matInfo.setVector(this->getStrain());
	case 3:
	{
		Vector res(3);
		res(0) = dmgtP;
		res(1) = dmgcP;
		res(2) = dmgsP;
		return matInfo.setVector(res);
	}
	case 4:
	{
		Vector res(4);
		res(0) = ktP;
		res(1) = kcP;
		res(2) = ktsP;
		res(3) = kcsP;
		return matInfo.setVector(res);
	}
	default:
		return -1;
	}
}

///////////  the methods for dealing geometry

const Vector
ASIConcrete::getVerticalLine(Vector v1)
{
	Vector line(3);
	line(0) = 1;
	line(2) = -v1(0);
	return line;
}

const Vector
ASIConcrete::getLineByTwoPoints(Vector v1, Vector v2)
{
	Vector para(3);

	if (v1 == v2)
		para(1) = 1.0;

	para(0) = (v2(1) - v1(1));
	para(1) = -(v2(0) - v1(0));
	para(2) = v1(1) * v2(0) - v1(0) * v2(1);

	if (para(1) < 0)
		para = -1 * para;

	return para;
}


void
ASIConcrete::getEllipseLineCrossPoint(Vector& crossP1, Vector& crossP2, Vector ellipse, Vector line)
{

	double const1, factor1, const2, factor2, delta, denominator;

	const1 = pow(ellipse(1) * line(1), 2) * ellipse(2) - ellipse(0) * ellipse(0) * line(0) * line(2);
	const2 = -ellipse(1) * ellipse(1) * line(1) * (ellipse(2) * line(0) + line(2));

	factor1 = -line(1);
	factor2 = line(0);

	delta = pow(ellipse(0) * ellipse(1), 2) * (pow(ellipse(0) * line(0), 2) - pow(ellipse(2) * line(0), 2) + pow(ellipse(1) * line(1), 2) -
		2 * ellipse(2) * line(0) * line(2) - line(2) * line(2));
	denominator = pow(ellipse(0) * line(0), 2) + pow(ellipse(1) * line(1), 2);

	if (delta < 0)
		delta = 0;

	crossP1(0) = (const1 + factor1 * pow(delta, 0.5)) / denominator;
	crossP1(1) = (const2 + factor2 * pow(delta, 0.5)) / denominator;
	crossP2(0) = (const1 - factor1 * pow(delta, 0.5)) / denominator;
	crossP2(1) = (const2 - factor2 * pow(delta, 0.5)) / denominator;

}


void
ASIConcrete::getTwoLinesCrossPoint(Vector& crossP, Vector line1, Vector line2)
{
	double numerator1, numerator2, denominator;

	numerator1 = line1(2) * line2(1) - line1(1) * line2(2);
	numerator2 = line1(0) * line2(2) - line1(2) * line2(0);
	denominator = line1(1) * line2(0) - line1(0) * line2(1);
	crossP(0) = numerator1 / denominator;
	crossP(1) = numerator2 / denominator;
}


////////////////////////  the mathods to determin dmgage value

void
ASIConcrete::getDmgTensile(const Vector line0, const Vector linei, Vector sig, double& dmgt)
{

	Vector direction(3), origin(2);
	Vector crossP0(2), crossPu(2);
	direction = getLineByTwoPoints(origin, sig);

	getTwoLinesCrossPoint(crossPu, direction, linei);

	if (sig.Norm() < crossPu.Norm())
	{
		getTwoLinesCrossPoint(crossP0, direction, line0);
		double dltDmg = (sig - crossP0).Norm() / (crossPu - crossP0).Norm();
		dmgt = dmgtP + (1 - dmgtP) * dltDmg;
	}
	else
	{
		dmgt = 1.0;
	}

	if (dmgt > 1.0)
		dmgt = 1.0;
}

void
ASIConcrete::getDmgShear(const Vector line0, const Vector linei, Vector sig, double& dmgt, double& dmgs)
{
	Vector direction(3), p1(2), p2(2);
	direction = getVerticalLine(sig);

	Vector crossP0(2), crossPi(2), crossInit(2);
	getTwoLinesCrossPoint(crossPi, direction, linei);

	if (abs(sig(1)) < abs(crossPi(1)))
	{
		if (eps(1) > 0)
		{
			getTwoLinesCrossPoint(crossP0, direction, line0);
		
			double dltDmg = (sig - crossP0).Norm() / (crossPi - crossP0).Norm();
			dmgt = dmgtP + (1 - dmgtP) * dltDmg * (1 - sig(0) / (fcyP / 2));
			dmgs = dmgsP + (1 - dmgsP) * dltDmg * sig(0) / (fcyP / 2);
		}
		else
		{
			getTwoLinesCrossPoint(crossP0, direction, line0);
		
			double dltDmg = (sig - crossP0).Norm() / (crossPi - crossP0).Norm();
			dmgt = dmgtP + (1 - dmgtP) * dltDmg * (1 - sig(0) / (fcyP / 2));
			dmgs = dmgsP + (1 - dmgsP) * dltDmg * sig(0) / (fcyP / 2);
		}
	}
	else
	{
		dmgt = dmgtP + (1 - dmgtP) * (1 - sig(0) / (fcyP / 2));
		dmgs = dmgsP + (1 - dmgsP) * sig(0) / (fcyP / 2);
	}

	if (dmgt > 1.0)
		dmgt = 1.0;

	if (dmgs > 1.0)
		dmgs = 1.0;

}

void
ASIConcrete::getDmgCompressive(const Vector ellipse0P, const Vector ellipseiP, Vector sig, double& dmgc, double& dmgs)
{
	Vector p(2), direction(3), ellipseInit(3);
	p(0) = ellipse0P(2);
	direction = getLineByTwoPoints(sig, p);

	Vector crossP01(2), crossP02(2), crossPi1(2), crossPi2(2), crossPInit1(2), crossPInit2(2);
	Vector crossP0(2), crossPi(2), crossPInit(2);

	getEllipseLineCrossPoint(crossPi1, crossPi2, ellipseiP, direction);
	if (((sig - p) ^ (crossPi1 - p)) > 0)
		crossPi = crossPi1;
	else
		crossPi = crossPi2;

	if (sig.Norm() < crossPi.Norm())
	{
		getEllipseLineCrossPoint(crossP01, crossP02, ellipse0P, direction);
		if (((sig - p) ^ (crossP01 - p)) > 0)
			crossP0 = crossP01;
		else
			crossP0 = crossP02;

		getEllipseLineCrossPoint(crossPInit1, crossPInit2, ellipseInit, direction);
		if (((sig - p) ^ (crossPInit1 - p)) > 0)
			crossPInit = crossPInit1;
		else
			crossPInit = crossPInit2;

		double dltDmg = (sig - crossP0).Norm() / (crossPi - crossP0).Norm();

		dmgc = dmgcP + (1 - dmgcP) * dltDmg * abs(p(0) - sig(0)) / (p - sig).Norm();
		dmgs = dmgsP + (1 - dmgsP) * dltDmg * abs(p(1) - sig(1)) / (p - sig).Norm();
	}
	else
	{
		dmgc = dmgcP + (1 - dmgcP) * abs(p(0) - sig(0)) / (p - sig).Norm();
		dmgs = dmgsP + (1 - dmgsP) * abs(p(1) - sig(1)) / (p - sig).Norm();
	}

	if (dmgc > ultmDmg)
		dmgc = ultmDmg;

	if (dmgs > 1.0)
		dmgs = 1.0;

}


////////////////////////  the mathods to calculate stress

void
ASIConcrete::getStressTensile(Vector sig, Vector& sigCor, double dmgt)
{
	Vector origin(2), direction(3);
	direction = getLineByTwoPoints(sig, origin);

	if (dmgt == 1.0)
		sigCor.Zero();
	else
	{
		Vector p1(2), p2(2), line(3);
		if (eps(1) > 0)
		{
			p1(0) = (1 - dmgt) * ft;
			p2(1) = (1 - dmgt) * fs;
			line = getLineByTwoPoints(p1, p2);
		}
		else
		{
			p1(0) = (1 - dmgt) * ft;
			p2(1) = -(1 - dmgt) * fs;
			line = getLineByTwoPoints(p1, p2);
		}
		 getTwoLinesCrossPoint(sigCor, line, direction);
	}
}

void
ASIConcrete::getStressShear(Vector sig, Vector& sigCor, double dmgc, double dmgts, double dmgcs)
{
	Vector direction(3);
	direction = getVerticalLine(sig);

	Vector p1(2), p2(2), line(3);

	double axRat = (aspectRatioUlti - aspectRatioInit) / ultmDmg * dmgcs + aspectRatioInit;

	if (sig(1) > 0)
	{
		p1(0) = (1 - dmgc) * fc / 2;
		p1(1) = abs(p1(0)) * axRat;
		p2(1) = (1 - dmgts) * fs;
		line = getLineByTwoPoints(p1, p2);
	}
	else
	{
		p1(0) = (1 - dmgc) * fc / 2;
		p1(1) = -abs(p1(0)) * axRat;
		p2(1) = -(1 - dmgts) * fs;
		line = getLineByTwoPoints(p1, p2);
	}
	getTwoLinesCrossPoint(sigCor, direction, line);
}


void
ASIConcrete::getStressCompressive(Vector sig, Vector& sigCor, double dmgc, double dmgcs)
{
	Vector p(2), direction(3);
	p(0) = (1 - dmgc) * fc / 2;
	direction = getLineByTwoPoints(sig, p);

	double axRat = (aspectRatioUlti - aspectRatioInit) / ultmDmg * dmgcs + aspectRatioInit;

	Vector ellipse(3);
	ellipse(2) = p(0);
	ellipse(0) = abs(ellipse(2));
	ellipse(1) = ellipse(0) * axRat;

	Vector sigTr1(2), sigTr2(2);
	getEllipseLineCrossPoint(sigTr1, sigTr2, ellipse, direction);
	if (((sig - p) ^ (sigTr1 - p)) > 0)
		sigCor = sigTr1;
	else
		sigCor = sigTr2;
}


///////////////////////////////  the mothods to calculate tangent

Matrix
ASIConcrete::getTangentTensile(Vector sig, Vector sigCor, const Vector line0P, const Vector lineiP)
{
	Vector epsTrN(2), epsTrS(2);
	epsTrN = eps;
	epsTrN(0) += perturb;

	epsTrS = eps;
	if (sig(1) > 0)
		epsTrS(1) += perturb;
	else
		epsTrS(1) -= perturb;

	Matrix tan(2, 2);
	Vector sigTrN(2), sigTrS(2);
	getTrialStressTangent(sigTrN, tan, epsTrN);
	getTrialStressTangent(sigTrS, tan, epsTrS);

	double dmgtPtb, dmgsPtb;
	Vector sigNTrCor(2), sigSTrCor(2);

	getDmgTensile(line0P, lineiP, sigTrN, dmgtPtb);
	getStressTensile(sigTrN, sigNTrCor, dmgtPtb);
	tan(0, 0) = (sigNTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = (sigNTrCor(1) - sigCor(1)) / perturb;

	getDmgTensile(line0P, lineiP, sigTrS, dmgtPtb);
	getStressTensile(sigTrS, sigSTrCor, dmgtPtb);
	if (sig(1) > 0)
	{
		tan(0, 1) = (sigSTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = (sigSTrCor(1) - sigCor(1)) / perturb;
	}
	else
	{
		tan(0, 1) = -(sigSTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = -(sigSTrCor(1) - sigCor(1)) / perturb;
	}

	if (abs(sig(1)) < 1e-4)
	{
		tan(1, 0) = 0.0;
		tan(0, 1) = 0.0;
	}

	if (abs(tan(0, 0)) < 1e-10)
		tan(0, 0) = 1.e-10;
	if (abs(tan(1, 1)) < 1e-10)
		tan(1, 1) = 1.e-10;

	return tan;
}

Matrix
ASIConcrete::getTangentShear(Vector sig, Vector sigCor, const Vector line0P, const Vector lineiP)
{
	Vector epsTrN(2), epsTrS(2);
	epsTrN = eps;
	epsTrN(0) += perturb;

	epsTrS = eps;
	if (sig(1) > 0)
		epsTrS(1) += perturb;
	else
		epsTrS(1) -= perturb;

	Matrix tan(2, 2);
	Vector sigTrN(2), sigTrS(2);
	getTrialStressTangent(sigTrN, tan, epsTrN);
	getTrialStressTangent(sigTrS, tan, epsTrS);

	double dmgsPtb, dmgtPtb;
	Vector sigTrCor(2);

	if (sigTrN(0) > 0)
	{
		getDmgTensile(line0P, lineiP, sigTrN, dmgtPtb);
		getStressTensile(sigTrN, sigTrCor, dmgtPtb);
	}
	else
	{
		getDmgShear(line0P, lineiP, sigTrS, dmgtPtb, dmgsPtb);
		getStressShear(sigTrS, sigTrCor, dmgcP, dmgtPtb, dmgsPtb);
	}
	tan(0, 0) = (sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = (sigTrCor(1) - sigCor(1)) / perturb;

	getDmgShear(line0P, lineiP, sigTrS, dmgtPtb, dmgsPtb);
	getStressShear(sigTrS, sigTrCor, dmgcP, dmgtPtb, dmgsPtb);
	if (sig(1) > 0)
	{
		tan(0, 1) = (sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = (sigTrCor(1) - sigCor(1)) / perturb;
	}
	else
	{
		tan(0, 1) = -(sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = -(sigTrCor(1) - sigCor(1)) / perturb;
	}

	if (abs(tan(0, 0)) < 1e-10)
		tan(0, 0) = 1.e-10;
	if (abs(tan(1, 1)) < 1e-10)
		tan(1, 1) = 1.e-10;

	return tan;
}


Matrix
ASIConcrete::getTangentShear(Vector sig, Vector sigCor, const Vector line0P)
{
	Vector epsTrN(2), epsTrS(2);
	epsTrN = eps;
	epsTrN(0) += perturb;

	epsTrS = eps;
	if (sig(1) > 0)
		epsTrS(1) += perturb;
	else
		epsTrS(1) -= perturb;

	Matrix tan(2, 2);
	Vector sigTrN(2), sigTrS(2);
	getTrialStressTangent(sigTrN, tan, epsTrN);
	getTrialStressTangent(sigTrS, tan, epsTrS);

	Vector sigTrCor(2);

	if (sigTrN(0) > 0)
		sigTrCor.Zero();
	else
		getStressShear(sigTrN, sigTrCor, dmgcP, 1.0, ultmDmg);
	tan(0, 0) = (sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = (sigTrCor(1) - sigCor(1)) / perturb;

	getStressShear(sigTrS, sigTrCor, dmgcP, 1.0, ultmDmg);
	if (sig(1) > 0)
	{
		tan(0, 1) = (sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = (sigTrCor(1) - sigCor(1)) / perturb;
	}
	else
	{
		tan(0, 1) = -(sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = -(sigTrCor(1) - sigCor(1)) / perturb;
	}

	if (abs(tan(0, 0)) < 1e-10)
		tan(0, 0) = 1.e-10;
	if (abs(tan(1, 1)) < 1e-10)
		tan(1, 1) = 1.e-10;

	return tan;
}



Matrix
ASIConcrete::getTangentCopressive(Vector sig, Vector sigCor, const Vector ellipse0P, const Vector ellipseiP)
{
	Vector epsTrN(2), epsTrS(2);
	epsTrN = eps;
	epsTrN(0) -= perturb;

	epsTrS = eps;
	if (sig(1) > 0)
		epsTrS(1) += perturb;
	else
		epsTrS(1) -= perturb;

	Matrix tan(2, 2);
	Vector sigTrN(2), sigTrS(2);
	getTrialStressTangent(sigTrN, tan, epsTrN);
	getTrialStressTangent(sigTrS, tan, epsTrS);

	double dmgcPtb, dmgsPtb;
	Vector sigTrCor(2);

	getDmgCompressive(ellipse0P, ellipseiP, sigTrN, dmgcPtb, dmgsPtb);
	getStressCompressive(sigTrN, sigTrCor, dmgcPtb, dmgsPtb);
	tan(0, 0) = -(sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = -(sigTrCor(1) - sigCor(1)) / perturb;

	getDmgCompressive(ellipse0P, ellipseiP, sigTrS, dmgcPtb, dmgsPtb);
	getStressCompressive(sigTrS, sigTrCor, dmgcPtb, dmgsPtb);
	if (sig(1) > 0)
	{
		tan(0, 1) = (sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = (sigTrCor(1) - sigCor(1)) / perturb;
	}
	else
	{
		tan(0, 1) = -(sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = -(sigTrCor(1) - sigCor(1)) / perturb;
	}

	if (abs(sig(1) < 1e-3))
		tan(0, 1) = 0;

	if (abs(tan(0, 0)) < 1e-10)
		tan(0, 0) = 1.e-10;
	if (abs(tan(1, 1)) < 1e-10)
		tan(1, 1) = 1.e-10;

	return tan;
}


Matrix
ASIConcrete::getTangentCopressive(Vector sig, Vector sigCor, const Vector ellipse0P)
{
	Vector epsTrN(2), epsTrS(2);
	epsTrN = eps;
	epsTrN(0) -= perturb;

	epsTrS = eps;
	if (sig(1) > 0)
		epsTrS(1) += perturb;
	else
		epsTrS(1) -= perturb;

	Matrix tan(2, 2);
	Vector sigTrN(2), sigTrS(2);
	getTrialStressTangent(sigTrN, tan, epsTrN);
	getTrialStressTangent(sigTrS, tan, epsTrS);

	Vector sigTrCor(2);
	getStressCompressive(sigTrN, sigTrCor, ultmDmg, ultmDmg);
	tan(0, 0) = -(sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = -(sigTrCor(1) - sigCor(1)) / perturb;

	getStressCompressive(sigTrS, sigTrCor, ultmDmg, ultmDmg);
	if (sig(1) > 0)
	{
		tan(0, 1) = (sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = (sigTrCor(1) - sigCor(1)) / perturb;
	}
	else
	{
		tan(0, 1) = -(sigTrCor(0) - sigCor(0)) / perturb;
		tan(1, 1) = -(sigTrCor(1) - sigCor(1)) / perturb;
	}

	if (abs(sig(1)) < 1e-4)
	{
		tan(1, 0) = 0.0;
		tan(0, 1) = 0.0;
	}

	if (abs(tan(0, 0)) < 1e-10)
		tan(0, 0) = 1.e-10;
	if (abs(tan(1, 1)) < 1e-10)
		tan(1, 1) = 1.e-10;

	return tan;
}



void
ASIConcrete::getTrialStressTangent(Vector& sigTr, Matrix& tangent, Vector epsTr)
{
	sigTr = sigP + eP * (epsTr - epsP);
	tangent = eP;

	if (epsTr(0) > epsresn)
	{
		if (epsTr(0) < epsP(0))
		{
			sigTr(0) = (epsTr(0) - epsresn) * sigP(0) / (epsP(0) - epsresn);
			tangent(0,0) = sigP(0) / (epsP(0) - epsresn);
		}
	}
	else if (epsTr(0) <= epsresn)
	{
		if (epsTr(0) > epsP(0))
		{
			sigTr(0) = (epsTr(0) - epsresn) * sigP(0) / (epsP(0) - epsresn);
			tangent(0, 0) = sigP(0) / (epsP(0) - epsresn);
		}
	}

	if (epsTr(1) < epsP(1) && epsTr(1) > 0)
	{
		sigTr(1) = sigP(1) / epsP(1) * epsTr(1);
		tangent(1, 1) = sigP(1) / epsP(1);
	}
	else if (epsTr(1) > epsP(1) && epsTr(1) < 0)
	{
		sigTr(1) = sigP(1) / epsP(1) * epsTr(1);
		tangent(1, 1) = sigP(1) / epsP(1);
	}
}

Vector 
ASIConcrete::getDamage()
{
	Vector dmg(3);
	dmg(0) = dmgtP;
	dmg(1) = dmgcP;
	dmg(2) = dmgsP;
	return dmg;
}
