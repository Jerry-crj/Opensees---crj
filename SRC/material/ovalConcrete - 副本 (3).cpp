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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ovalConcrete.cpp,v $

// Written: fmk
// Created: 03/06
//
// Description: This file contains the class implementation of ovalConcrete. 
// This ovalConcrete is based on an f2c of the FEDEAS material
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

#include <ovalConcrete.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <MaterialResponse.h>

#define MAT_TAG_OVALConcrete 123123141232

void*
OPS_ovalConcrete()
{
	// Pointer to a uniaxial material that will be returned
	NDMaterial* theMaterial = 0;

	int    iData[1];
	double dData[13];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid NDMaterial ovalConcrete tag" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData < 10) {
		opserr << "Invalid #args, want: NDMaterial ovalConcrete " << iData[0] << " fpc? epsc0? epscu? epst0? epstu?\n";
		return 0;
	}

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid #args, want: NDMaterial ovalConcrete " << iData[0] << " fpc? epsc0? epscu? epst0? epstu?\n";
		return 0;
	}

	// Parsing was successful, allocate the material
	if (numData == 10)
	{
		theMaterial = new ovalConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9]);
	}
	else if (numData == 11)
	{
		theMaterial = new ovalConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9], dData[10]);
	}
	else if (numData == 12)
	{
		theMaterial = new ovalConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);
	}
	else if (numData == 13)
	{
		theMaterial = new ovalConcrete(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
			dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]);
	}

	if (theMaterial == 0) {
		opserr << "WARNING could not create NDMaterial of type ovalConcrete Material\n";
		return 0;
	}

	return theMaterial;
}

ovalConcrete::ovalConcrete(int tag, double sc0, double ec0, double scu, double ecu, double rat,
	double st, double ets, double ss0, double grat, double gs, double rmax, double rmin, double p) :
	NDMaterial(tag, MAT_TAG_OVALConcrete), dmgt(0), dmgtP(0), dmgc(0), dmgcP(0),
	dmgcsp(0), dmgcspP(0), dmgcsn(0), dmgcsnP(0), dmgtsp(0), dmgtspP(0), dmgtsn(0), dmgtsnP(0),
	epsresn(0), epsresnMin(0), epsressp(0), epsresspMin(0), epsressn(0), epsressnMin(0),
	eps(2), epsP(2), sig(2), sigCor(2), sigP(2), e(2, 2), e0(2, 2), eP(2, 2), unldrat(rat), rho(0)
{
	Ets = ets;
	Gs = gs;

	ft0init = abs(st);
	fs0init = abs(ss0);
	fc0init = -abs(sc0);
	fcuinit = -abs(scu);

	kcP = fc0init / ec0;
	ktP = kcP;
	ktspP = kcP * grat;
	ktsnP = kcP * grat;
	kcspP = ktspP;
	kcsnP = ktsnP;
	knInit = kcP;
	ksInit = ktspP;

	epsc0 = -abs(ec0);
	epscu = -abs(ecu);
	epst0 = abs(st / kcP);
	epstu = epst0 + abs(st / Ets);
	epss0 = abs(ss0 / ktspP);
	epssu = epss0 + abs(ss0 / gs);

	ft0 = ft0init;
	fc0 = fc0init;
	fs0 = fs0init;

	ft0P = abs(ft0);
	fsp0P = abs(fs0);
	fsn0P = abs(fs0);
	fc0P = -abs(fc0);
	ftiP = ft0P + ktP * (epstu - epst0);
	fciP = fc0P + kcP * (epscu - epsc0);
	ftspiP = fsp0P + ktspP * (epssu - epss0);
	ftsniP = fsn0P + ktsnP * (epssu - epss0);
	fcspiP = fsp0P + kcspP * (epssu - epss0);
	fcsniP = fsn0P + kcsnP * (epssu - epss0);

	ultmDmg = 1 - fcuinit / fc0init;

	eP(0, 0) = ktP;
	eP(1, 1) = ktspP;

	epsp = (ecu * kcP * unldrat - scu) / (kcP * (-1 + unldrat));
	sigp = (ecu * kcP * unldrat - scu) / (-1 + unldrat);

	perturb = p;

	axRatpP = rmax;
	axRatnP = rmax;
	axRatMax = rmax;
	axRatMin = rmin;
}

ovalConcrete::ovalConcrete(void) :
	NDMaterial(0, MAT_TAG_OVALConcrete),
	fc0(0), epsc0(0), epscu(0), epst0(0), epstu(0), epss0(0), epssu(0),
	fc0init(0), ft0init(0), fs0init(0), fcuinit(0), ft0(0), fs0(0),
	fc0P(0), ft0P(0), fsp0P(0), fsn0P(0), fciP(0), ftiP(0), ftspiP(0), ftsniP(0), fcspiP(0), fcsniP(0),
	dmgt(0), dmgtP(0), dmgc(0), dmgcP(0), dmgcsp(0), dmgcspP(0), dmgcsn(0), dmgcsnP(0),
	dmgtsp(0), dmgtspP(0), dmgtsn(0), dmgtsnP(0),
	ktspP(0), ktsnP(0), kcspP(0), kcsnP(0), kcP(0), ktP(0), ksInit(0), knInit(0),
	epsresn(0), epsresnMin(0), epsressp(0), epsresspMin(0), epsressn(0), epsressnMin(0),
	Ets(0), Gs(0), epsP(2), sigP(2), eps(2), sig(2), e(2, 2), eP(2, 2),
	perturb(1e-7), axRatpP(0.3), axRatnP(0.3), axRatMax(0.3), axRatMin(0.3),
	epsp(0.0), sigp(0.0), unldrat(0), ultmDmg(1.0), rho(0)
{

}

ovalConcrete::~ovalConcrete(void)
{
	// Does nothing
}

NDMaterial*
ovalConcrete::getCopy(void)
{
	ovalConcrete* theCopy = new ovalConcrete(this->getTag(), fc0init, epsc0, fcuinit, epscu,
		unldrat, ft0init, Ets, fs0init, ktspP / ktP, Gs, axRatMax, axRatMin, perturb);

	return theCopy;
}

NDMaterial*
ovalConcrete::getCopy(const char* type)
{
	ovalConcrete* theCopy = 0;
	theCopy = (ovalConcrete*)this->getCopy();

	return theCopy;
}


int
ovalConcrete::setTrialStrain(const Vector& trialStrain)
{
	eps = trialStrain;

	getTrialStressTangent(sig, e, eps);

	if (dmgtP > 1.0 - 1.0e-10 && sig(0) > 0)
	{
		dmgt = 1.0;
		sig(0) = 0;
	}

	if (sig(0) > 0)
	{
		if (dmgtspP > 1.0 - 1.0e-10 && sig(1) > 0)
		{
			dmgtsp = 1.0;
			sig(1) = 0;
		}if (dmgtsnP > 1.0 - 1.0e-10 && sig(1) < 0)
		{
			dmgtsn = 1.0;
			sig(1) = 0;
		}
	}
	else
	{
		if (dmgcspP > 1.0 - 1.0e-10 && sig(1) > 0)
		{
			dmgcsp = 1.0;
			sig(1) = 0;
		}if (dmgcsnP > 1.0 - 1.0e-10 && sig(1) < 0)
		{
			dmgcsn = 1.0;
			sig(1) = 0;
		}
	}

	Vector sigCor(2);

	if (sig(0) > 0)
	{
		if (dmgtP == 1.0 && dmgtspP == 1.0 && dmgtsnP == 1.0)
		{
			sig.Zero();
			e.Zero();
			e(0, 0) = 1.0e-10;
			e(1, 1) = 1.0e-10;
		}
		if (sig(1) >= 0)
		{
			if (dmgtP == 1.0 && dmgtspP == 1.0)
			{
				sig.Zero();
				e.Zero();
				e(0, 0) = 1.0e-10;
				e(1, 1) = 1.0e-10;
			}
			else if (dmgtP == 1.0 && dmgtspP < 1.0)
			{
				sig(0) = 0.0;
				e(0, 0) = 1.0e-10;
				if (sig(1) < fsn0P)
				{
					dmgtsp = dmgtspP + (1 - dmgtspP) * (sig(1) - fsn0P) / (ftsniP - fsn0P);

					sig(1) = fs0init * (1 - dmgtsp);

					double dmgPtb = dmgtspP + (1 - dmgtspP) * (sig(1) + perturb - fsn0P) / (ftsniP - fsn0P);

					double sigPtb = fs0init * (1 - dmgPtb);

					e(1, 1) = (sig(1) - sigPtb) / perturb;
				}
			}
			else if (dmgtP < 1.0 && dmgtspP == 1.0)
			{
				sig(1) = 0.0;
				e(1, 1) = 1.0e-10;

				dmgt = dmgtP + (1 - dmgtP) * (sig(0) - ft0P) / (ftiP - ft0P);

				sig(0) = ft0init * (1 - dmgt);

				double dmgPtb = dmgtP + (1 - dmgtP) * (sig(1) + perturb - ft0P) / (ftiP - ft0P);

				double sigPtb = ft0init * (1 - dmgPtb);

				e(1, 1) = (sig(0) - sigPtb) / perturb;
			}
			else
			{
				if (sig(0) / ft0P + abs(sig(1) / fsp0P) - 1 > 0.0)
				{
					if (sig(0) / ftiP + abs(sig(1) / ftspiP) - 1 < 0.0)
					{
						Vector p1(2), p2(2), line0P(3), lineiP(3);

						p1(0) = ft0P;
						p2(1) = fsp0P;
						line0P = getLineByTwoPoints(p1, p2);
						p1(0) = ftiP;
						p2(1) = ftspiP;
						lineiP = getLineByTwoPoints(p1, p2);

						double dDmgt, dDmgs;
						getDltDmgTensile(line0P, lineiP, sig, dDmgt, dDmgs);

						dmgt = dmgtP + (1 - dmgtP) * dDmgt;
						dmgtsp = dmgtspP + (1 - dmgtspP) * dDmgs;
						double minDmgt = min(dmgtsp, dmgtsnP);
						if (dmgt < minDmgt)
							dmgt = minDmgt;
						if (dmgt > 1.0)
							dmgt = 1.0;
						if (dmgtsp > 1.0)
							dmgtsp = 1.0;

						getStressTensile(sig, sigCor, dmgt, dmgtsp);

						e = getTangentTensile(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
					else
					{
						Vector p1(2), p2(2), line0P(3), lineiP(3);

						p1(0) = ft0P;
						p2(1) = fsp0P;
						line0P = getLineByTwoPoints(p1, p2);
						p1(0) = ftiP;
						p2(1) = ftspiP;
						lineiP = getLineByTwoPoints(p1, p2);

						double dDmg = min(1 - dmgtP, 1 - dmgtspP);

						dmgt = dmgtP + dDmg * sig(0) / sig.Norm();
						dmgtsp = dmgtspP + dDmg * abs(sig(1)) / sig.Norm();

						double minDmgt = min(dmgtsp, dmgtsnP);
						if (dmgt < minDmgt)
							dmgt = minDmgt;
						if (dmgt > 1.0)
							dmgt = 1.0;
						if (dmgtsp > 1.0)
							dmgtsp = 1.0;

						getStressTensile(sig, sigCor, dmgt, dmgtsp);

						e = getTangentTensile(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
				}
			}
		}
		else if (sig(1) <= 0)
		{
			if (dmgtP == 1.0 && dmgtsnP == 1.0)
			{
				sig.Zero();
				e.Zero();
				e(0, 0) = 1.0e-10;
				e(1, 1) = 1.0e-10;
			}
			else if (dmgtP == 1.0 && dmgtsnP < 1.0)
			{
				sig(0) = 0.0;
				e(0, 0) = 1.0e-10;
				if (sig(1) > fsn0P)
				{
					dmgtsn = dmgtsnP + (1 - dmgtsnP) * (sig(1) - fsn0P) / (ftsniP - fsn0P);

					sig(1) = fs0init * (1 - dmgtsn);

					double dmgPtb = dmgtsnP + (1 - dmgtsnP) * (sig(1) + perturb - fsn0P) / (ftsniP - fsn0P);

					double sigPtb = fs0init * (1 - dmgPtb);

					e(1, 1) = (sig(1) - sigPtb) / perturb;
				}
			}
			else if (dmgtP < 1.0 && dmgtsnP == 1.0)
			{
				sig(1) = 0.0;
				e(1, 1) = 1.0e-10;

				dmgt = dmgtP + (1 - dmgtP) * (sig(0) - ft0P) / (ftiP - ft0P);

				sig(0) = ft0init * (1 - dmgt);

				double dmgPtb = dmgtP + (1 - dmgtP) * (sig(1) + perturb - ft0P) / (ftiP - ft0P);

				double sigPtb = ft0init * (1 - dmgPtb);

				e(1, 1) = (sig(0) - sigPtb) / perturb;
			}
			else
			{
				if (sig(0) / ft0P + abs(sig(1) / fsn0P) - 1 > 0.0)
				{
					if (sig(0) / ftiP + abs(sig(1) / ftsniP) - 1 < 0.0)
					{
						Vector p1(2), p2(2), line0P(3), lineiP(3);

						p1(0) = ft0P;
						p2(1) = -fsn0P;
						line0P = getLineByTwoPoints(p1, p2);
						p1(0) = ftiP;
						p2(1) = -ftsniP;
						lineiP = getLineByTwoPoints(p1, p2);

						double dDmgt, dDmgs;
						getDltDmgTensile(line0P, lineiP, sig, dDmgt, dDmgs);

						dmgt = dmgtP + (1 - dmgtP) * dDmgt;
						dmgtsn = dmgtsnP + (1 - dmgtsnP) * dDmgs;

						double minDmgt = min(dmgtsn, dmgtspP);
						if (dmgt < minDmgt)
							dmgt = minDmgt;
						if (dmgt > 1.0)
							dmgt = 1.0;
						if (dmgtsn > 1.0)
							dmgtsn = 1.0;

						getStressTensile(sig, sigCor, dmgt, dmgtsn);

						e = getTangentTensile(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
					else
					{
						Vector p1(2), p2(2), line0P(3), lineiP(3);

						p1(0) = ft0P;
						p2(1) = -fsn0P;
						line0P = getLineByTwoPoints(p1, p2);
						p1(0) = ftiP;
						p2(1) = -ftsniP;
						lineiP = getLineByTwoPoints(p1, p2);

						double dDmg = min(1 - dmgtP, 1 - dmgtsnP);

						dmgt = dmgtP + dDmg * sig(0) / sig.Norm();
						dmgtsn = dmgtsnP + dDmg * abs(sig(1)) / sig.Norm();

						double minDmgt = min(dmgtsn, dmgtspP);
						if (dmgt < minDmgt)
							dmgt = minDmgt;
						if (dmgt > 1.0)
							dmgt = 1.0;
						if (dmgtsn > 1.0)
							dmgtsn = 1.0;

						getStressTensile(sig, sigCor, dmgt, dmgtsn);

						e = getTangentTensile(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
				}
			}
		}
	}
	else if (sig(0) > fc0P / 2)
	{
		Vector p1(2), p2(2), line0P(3), lineiP(3);

		if (sig(1) > 0)  // 正值
		{
			p1(0) = fc0P / 2;
			p1(1) = abs(fc0P) / 2 * axRatpP;
			p2(1) = fsp0P;
			line0P = getLineByTwoPoints(p1, p2);
			if (line0P(0) * sig(0) + line0P(1) * sig(1) + line0P(2) > 0) // out of yield surface
			{
				if (dmgtspP > 1.0 - 1.0e-10 && dmgcspP > ultmDmg - 1.0e-10)
				{
					dmgtsp = 1.0;
					dmgcsp = ultmDmg;

					getStressShear(sig, sigCor, dmgcP, 1.0, ultmDmg);

					e = getTangentShear(sig, sigCor, line0P);

					sig = sigCor;
				}
				else
				{
					p1(0) = fc0P / 2;
					p1(1) = abs(fciP - fc0P / 2) * axRatpP;
					p2(1) = fcspiP;
					lineiP = getLineByTwoPoints(p1, p2);

					if (lineiP(0) * sig(0) + lineiP(1) * sig(1) + lineiP(2) < 0) // in ultimate surface
					{
						double dDmg = getDltDmgShear(line0P, lineiP, sig);

						dmgcsp = dmgcspP + (1 - dmgcspP) * dDmg * sig(0) / fc0P / 2;
						//dmgtsp = dmgtspP + (1 - dmgtspP) * dDmg;
						dmgtsp = dmgtspP + dDmg;

						double minDmgt = min(dmgtsp, dmgtsnP);
						if (dmgt < minDmgt)
							dmgt = minDmgt;
						if (dmgcsp > ultmDmg)
						{
							dmgcsp = ultmDmg;
							dmgtsp = dmgtspP + (1 - ultmDmg) * dDmg;
						}
						if (dmgtsp > 1.0)
							dmgtsp = 1.0;

						getStressShear(sig, sigCor, dmgcP, dmgtsp, dmgcsp);

						e = getTangentShear(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
					else
					{
						double dDmg;
						dDmg = min(1.0 - dmgtsp, ultmDmg - dmgcsp);
						dmgcsp = dmgcspP + dDmg * sig(0) / fc0P / 2;
						dmgtsp = dmgtspP + dDmg;

						getStressShear(sig, sigCor, dmgcP, dmgtsp, dmgcsp);

						e = getTangentShear(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
				}
			}
		}
		else // 负值
		{
			p1(0) = fc0P / 2;
			p1(1) = -abs(fc0P) / 2 * axRatnP;
			p2(1) = -fsn0P;
			line0P = getLineByTwoPoints(p1, p2);

			if (line0P(0) * sig(0) + line0P(1) * sig(1) + line0P(2) < 0) // out of yield surface
			{
				if (dmgtsnP > 1.0 - 1.0e-10 && dmgcsnP > ultmDmg - 1.0e-10)
				{
					dmgtsn = 1.0;
					dmgcsn = ultmDmg;

					getStressShear(sig, sigCor, dmgcP, 1.0, ultmDmg);

					e = getTangentShear(sig, sigCor, line0P);

					sig = sigCor;
				}
				else
				{
					p1(0) = fc0P / 2;
					p1(1) = -abs(fciP - fc0P / 2) * axRatnP;
					p2(1) = -fcsniP;
					lineiP = getLineByTwoPoints(p1, p2);

					if (lineiP(0) * sig(0) + lineiP(1) * sig(1) + lineiP(2) > 0) // in ultimate surface
					{
						double dDmg = getDltDmgShear(line0P, lineiP, sig);

						dmgcsn = dmgcsnP + (1 - dmgcsnP) * dDmg * sig(0) / fc0P / 2;
						//dmgtsn = dmgtsnP + (1 - dmgtsnP) * dDmg;
						dmgtsn = dmgtsnP + dDmg;

						double minDmgt = min(dmgtsn, dmgtspP);
						if (dmgt < minDmgt)
							dmgt = minDmgt;
						if (dmgcsn > ultmDmg)
						{
							dmgcsn = ultmDmg;
							dmgtsn = dmgtsnP + (1 - ultmDmg) * dDmg;
						}
						if (dmgtsn > 1.0)
							dmgtsn = 1.0;

						getStressShear(sig, sigCor, dmgcP, dmgtsn, dmgcsn);

						e = getTangentShear(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
					else  // bug here
					{
						double dDmg;
						dDmg = min(1.0 - dmgtsn, ultmDmg - dmgcsn);
						dmgcsn = dmgcsnP + dDmg * sig(0) / fc0P / 2;
						dmgtsn = dmgtsnP + dDmg;

						getStressShear(sig, sigCor, dmgcP, dmgtsn, dmgcsn);

						e = getTangentShear(sig, sigCor, line0P, lineiP);

						sig = sigCor;
					}
				}
			}
		}
	}
	else   // shear and compressive area
	{
		if (sig(1) > 0)
		{
			Vector ellipse0P(3);
			ellipse0P(0) = abs(fc0P) / 2;
			ellipse0P(1) = ellipse0P(0) * axRatpP;
			ellipse0P(2) = fc0P / 2;

			if (pow((sig(0) - ellipse0P(2)) / ellipse0P(0), 2) + pow(sig(1) / ellipse0P(1), 2) - 1 > 0.0)  // outside of the ultimate sufrace
			{
				if (dmgcP > ultmDmg - 1.0e-10)
				{
					getStressCompressive(sig, sigCor, ultmDmg, ultmDmg);

					e = getTangentCopressive(sig, sigCor, ellipse0P);

					sig = sigCor;
				}
				else
				{
					Vector ellipseiP(3);
					ellipseiP(0) = abs(fciP - fc0P / 2);
					ellipseiP(1) = ellipseiP(0) * axRatpP;
					ellipseiP(2) = fc0P / 2;

					if (pow((sig(0) - ellipseiP(2)) / ellipseiP(0), 2) + pow(sig(1) / ellipseiP(1), 2) - 1 < 0.0)
					{
						double dDmg = getDltDmgCompressive(ellipse0P, ellipseiP, sig);

						dmgc = dmgcP + (1 - dmgcP) * dDmg;
						dmgt = dmgtP + (1 - dmgtP) * dDmg;
						dmgcsp = dmgcspP + (1 - dmgcspP) * dDmg;
						dmgcsn = dmgcsnP + (1 - dmgcsnP) * dDmg;

						if (dmgc > ultmDmg)
							dmgc = ultmDmg;
						if (dmgt > 1.0)
							dmgt = 1.0;
						if (dmgcsp > ultmDmg)
							dmgcsp = ultmDmg;
						if (dmgcsn > ultmDmg)
							dmgcsn = ultmDmg;

						if (sig(1) > 0)
							getStressCompressive(sig, sigCor, dmgc, dmgcsp);
						else
							getStressCompressive(sig, sigCor, dmgc, dmgcsp);

						e = getTangentCopressive(sig, sigCor, ellipse0P, ellipseiP);

						sig = sigCor;
					}
					else
					{
						dmgc = ultmDmg;
						dmgcsp = ultmDmg;
						dmgcsn = ultmDmg;

						getStressCompressive(sig, sigCor, ultmDmg, ultmDmg);

						e = getTangentCopressive(sig, sigCor, ellipse0P, ellipseiP);

						sig = sigCor;
					}
				}
			}
		}
		else
		{
			Vector ellipse0P(3);
			ellipse0P(0) = abs(fc0P) / 2;
			ellipse0P(1) = ellipse0P(0) * axRatnP;
			ellipse0P(2) = fc0P / 2;

			if (pow((sig(0) - ellipse0P(2)) / ellipse0P(0), 2) + pow(sig(1) / ellipse0P(1), 2) - 1 > 0.0)  // outside of the ultimate sufrace
			{
				if (dmgcP > ultmDmg - 1.0e-10)
				{
					getStressCompressive(sig, sigCor, ultmDmg, ultmDmg);

					e = getTangentCopressive(sig, sigCor, ellipse0P);

					sig = sigCor;
				}
				else
				{
					Vector ellipseiP(3);
					ellipseiP(0) = abs(fciP - fc0P / 2);
					ellipseiP(1) = ellipseiP(0) * axRatnP;
					ellipseiP(2) = fc0P / 2;

					if (pow((sig(0) - ellipseiP(2)) / ellipseiP(0), 2) + pow(sig(1) / ellipseiP(1), 2) - 1 < 0.0)
					{
						double dDmg = getDltDmgCompressive(ellipse0P, ellipseiP, sig);

						dmgc = dmgcP + (1 - dmgcP) * dDmg;
						dmgt = dmgtP + (1 - dmgtP) * dDmg;
						dmgcsp = dmgcspP + (1 - dmgcspP) * dDmg;
						dmgcsn = dmgcsnP + (1 - dmgcsnP) * dDmg;

						if (dmgc > ultmDmg)
							dmgc = ultmDmg;
						if (dmgt > 1.0)
							dmgt = 1.0;
						if (dmgcsp > ultmDmg)
							dmgcsp = ultmDmg;
						if (dmgcsn > ultmDmg)
							dmgcsn = ultmDmg;

						if (sig(1) > 0)
							getStressCompressive(sig, sigCor, dmgc, dmgcsp);
						else
							getStressCompressive(sig, sigCor, dmgc, dmgcsn);

						e = getTangentCopressive(sig, sigCor, ellipse0P, ellipseiP);

						sig = sigCor;
					}
					else
					{
						dmgc = ultmDmg;
						dmgcsp = ultmDmg;
						dmgcsn = ultmDmg;

						getStressCompressive(sig, sigCor, ultmDmg, ultmDmg);

						e = getTangentCopressive(sig, sigCor, ellipse0P, ellipseiP);

						sig = sigCor;
					}
				}
			}
		}
	}
	if (dmgtP > 1.0)
		bool debug = true;
	if (dmgcsp > 0 || dmgcsn > 0)
		bool debug = true;



	return 0;
}


int
ovalConcrete::setTrialStrain(const Vector& v, const Vector& r)
{
	return setTrialStrain(v);
}

int
ovalConcrete::setTrialStrainIncr(const Vector& v)
{
	Vector temp(3);
	temp = eps + v;
	return setTrialStrain(temp);
}

int
ovalConcrete::setTrialStrainIncr(const Vector& v, const Vector& rate)
{
	Vector temp(3);
	temp = eps + v;
	return setTrialStrain(temp);
}

const Vector&
ovalConcrete::getStrain(void)
{
	return eps;
}

const Vector&
ovalConcrete::getStress(void)
{
	return sig;
}

const Matrix&
ovalConcrete::getTangent(void)
{
	return e;
}

double
ovalConcrete::getRho(void)
{
	return rho;
}

const Matrix&
ovalConcrete::getInitialTangent(void)
{
	return e0;
}


int
ovalConcrete::commitState(void)
{
	sigP = sig;
	epsP = eps;

	ft0P = ft0init * (1 - dmgt);
	fc0P = fc0init * (1 - dmgc);
	fsp0P = fs0init * (1 - dmgtsp);
	fsn0P = fs0init * (1 - dmgtsn);

	ktP = ft0init * (1 - dmgt) / (epst0 + dmgt * (epstu - epst0));
	kcP = (fc0init * (1 - dmgc) - sigp) / (epsc0 + dmgc * (epscu - epsc0) / ultmDmg - epsp);
	ktspP = fs0init * (1 - dmgtsp) / (epss0 + dmgtsp * (epssu - epss0));
	ktsnP = fs0init * (1 - dmgtsn) / (epss0 + dmgtsn * (epssu - epss0));
	kcspP = fs0init * (1 - dmgcsp) / (epss0 + dmgcsp * (epssu - epss0));
	kcsnP = fs0init * (1 - dmgcsn) / (epss0 + dmgcsn * (epssu - epss0));

	ftiP = ft0P + ktP * (epstu - epst0) * (1 - dmgt);
	fciP = fc0P + kcP * (epscu - epsc0) / ultmDmg * (1 - dmgc);
	ftspiP = fsp0P + ktspP * (epssu - epss0) * (1 - dmgtsp);
	ftsniP = fsn0P + ktsnP * (epssu - epss0) * (1 - dmgtsn);
	fcspiP = fsp0P + kcspP * (epssu - epss0) * (1 - dmgcsp);
	fcsniP = fsn0P + kcsnP * (epssu - epss0) * (1 - dmgcsn);

	if (dmgc > dmgcP)
	{
		if (epsresnMin > eps(0) - sig(0) / kcP)
			epsresnMin = eps(0) - sig(0) / kcP;
	}
	else if (dmgc == ultmDmg)
	{
		if (sig(0) < fc0P / 2 &&
			pow((sig(0) - fc0P / 2) / abs(fc0P / 2), 2) + pow(sig(1) / abs(fc0P / 2) / axRatpP, 2) - 1 > -1.0e-10)
			if (epsresnMin > eps(0) - sig(0) / kcP)
				epsresnMin = eps(0) - sig(0) / kcP;
	}
	epsresn = epsresnMin;

	if (sig(1) > 0)
	{
		if (sig(0) > 0 && dmgtsp > dmgtspP && ktspP > 1.0e-3)
		{
			if (epsresspMin < eps(1) - sig(1) / ktspP)
				epsresspMin = eps(1) - sig(1) / ktspP;
		}
		else if (sig(0) <= 0)
		{
			if (dmgcsp > dmgcspP || dmgtsp > dmgtspP)
			{
				epsresspMin = eps(1) - sig(1) / kcspP;
			}
			else if (dmgcspP == ultmDmg && sig(1) > -axRatMin * sig(0))
			{
				if (epsresspMin < eps(1) - sig(1) / kcspP)
					epsresspMin = eps(1) - sig(1) / kcspP;
			}

			else if (dmgcspP == ultmDmg && pow((sig(0) - fc0P / 2) / abs(fc0P / 2), 2) + pow(sig(1) / abs(fc0P / 2) / axRatpP, 2) - 1 > -1.0e-10)
			{
				if (epsresspMin < eps(1) - sig(1) / kcspP)
					epsresspMin = eps(1) - sig(1) / kcspP;
			}
		}
	}
	else
	{
		if (sig(0) > 0 && dmgtsn > dmgtsnP && ktsnP > 1.0e-3)
		{
			if (epsressnMin > eps(1) - sig(1) / ktsnP)
				epsressnMin = eps(1) - sig(1) / ktsnP;
		}
		else if (sig(0) <= 0)
		{
			if (dmgcsn > dmgcsnP || dmgtsn > dmgtsnP)
			{
				if (epsressnMin > eps(1) - sig(1) / kcsnP)
					epsressnMin = eps(1) - sig(1) / kcsnP;
			}
			else if (dmgcsnP == ultmDmg && sig(1) < axRatMin * sig(0))
			{
				if (epsressnMin > eps(1) - sig(1) / kcsnP)
					epsressnMin = eps(1) - sig(1) / kcsnP;
			}
			else if (dmgcsnP == ultmDmg && pow((sig(0) - fc0P / 2) / abs(fc0P / 2), 2) + pow(sig(1) / abs(fc0P / 2) / axRatpP, 2) - 1 > -1.0e-10)
			{
				if (epsressnMin > eps(1) - sig(1) / kcsnP)
					epsressnMin = eps(1) - sig(1) / kcsnP;
			}
		}
	}
	epsressp = epsresspMin;
	epsressn = epsressnMin;

	if (eps(0) > epsresn)
	{
		eP(0, 0) = ktP;
		if (eps(1) > epsressp)
			eP(1, 1) = ktspP;
		else if (eps(1) < epsressn)
			eP(1, 1) = ktsnP;
		else
			eP(1, 1) = 1.0e-10;
	}
	else
	{
		eP(0, 0) = kcP;
		if (eps(1) > epsressp)
			eP(1, 1) = kcspP;
		else if (eps(1) < epsressn)
			eP(1, 1) = kcsnP;
		else
			eP(1, 1) = 1.0e-10;
	}

	dmgtP = dmgt;
	dmgcP = dmgc;
	dmgcspP = dmgcsp;
	dmgcsnP = dmgcsn;
	dmgtspP = dmgtsp;
	dmgtsnP = dmgtsn;

	axRatpP = axRatMax * (1 - dmgcsp) + axRatMin * dmgcsp;
	axRatnP = axRatMax * (1 - dmgcsn) + axRatMin * dmgcsn;

	return 0;
}

int
ovalConcrete::revertToLastCommit(void)
{
	dmgt = dmgtP;
	dmgc = dmgcP;
	dmgcsn = dmgcsnP;
	dmgcsp = dmgcspP;
	dmgtsn = dmgtsnP;
	dmgtsp = dmgtspP;

	sig = sigP;
	eps = epsP;

	return 0;
}

int
ovalConcrete::revertToStart(void)
{
	dmgt = 0.0;
	dmgc = 0.0;
	dmgcsp = 0.0;
	dmgcsn = 0.0;
	dmgtP = 0.0;
	dmgcP = 0.0;
	dmgcspP = 0.0;
	dmgcsnP = 0.0;

	ft0P = ft0init;
	fsp0P = fs0init;
	fsn0P = fs0init;
	fc0P = fc0init;
	ftiP = ft0P + ktP * (epstu - epst0);
	fciP = fc0P + kcP * (epscu - epsc0);
	ftspiP = fsp0P + ktspP * (epssu - epss0);
	ftsniP = fsn0P + ktsnP * (epssu - epss0);
	fcspiP = fsp0P + kcspP * (epssu - epss0);
	fcsniP = fsn0P + kcsnP * (epssu - epss0);

	sig.Zero();
	eps.Zero();
	sigP.Zero();
	epsP.Zero();

	return 0;
}

int
ovalConcrete::sendSelf(int commitTag, Channel& theChannel)
{
	return 0;
}

int
ovalConcrete::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
}

void
ovalConcrete::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "ovalConcrete:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ovalConcrete\", ";
		s << "\"Ec\": " << fc0 / epsc0 << ", ";
		s << "\"fc\": " << fc0 << ", ";
		s << "\"epsc\": " << epsc0 << ", ";
		s << "\"epscu\": " << epscu << ", ";
	}
}

Response*
ovalConcrete::setResponse(const char** argv, int argc,
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
ovalConcrete::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case 1:
		return matInfo.setVector(this->getStress());
	case 2:
		return matInfo.setVector(this->getStrain());
	case 3:
	{
		Vector res(6);
		res(0) = dmgtP;
		res(1) = dmgcP;
		res(2) = dmgtspP;
		res(3) = dmgcspP;
		res(4) = dmgtsnP;
		res(5) = dmgcsnP;
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

///////////  the methods for dealing geometry

const Vector
ovalConcrete::getVerticalLine(Vector v1)
{
	Vector line(3);
	line(0) = 1;
	line(2) = -v1(0);
	return line;
}

const Vector
ovalConcrete::getLineByTwoPoints(Vector v1, Vector v2)
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
ovalConcrete::getEllipseLineCrossPoint(Vector& crossP1, Vector& crossP2, Vector ellipse, Vector line)
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
ovalConcrete::getTwoLinesCrossPoint(Vector& crossP, Vector line1, Vector line2)
{
	double numerator1, numerator2, denominator;

	numerator1 = line1(2) * line2(1) - line1(1) * line2(2);
	numerator2 = line1(0) * line2(2) - line1(2) * line2(0);
	denominator = line1(1) * line2(0) - line1(0) * line2(1);
	crossP(0) = numerator1 / denominator;
	crossP(1) = numerator2 / denominator;
}

void
ovalConcrete::getPerpendicularPoint(Vector& objPoint, Vector line, Vector point)
{
	objPoint(0) = -(line(0) * line(2) - line(1) * line(1) * point(0) + line(0) * line(1) * point(1)) / (line(0) * line(0) + line(1) * line(1));
	objPoint(1) = -(line(1) * line(2) + line(0) * line(1) * point(0) - line(0) * line(0) * point(1)) / (line(0) * line(0) + line(1) * line(1));

	if (objPoint(1) * point(1) <0)
	{
		objPoint(0) = -line(2) / line(0);
		objPoint(1) = 0;
	}
	if (objPoint(0) < 0)
	{
		objPoint(0) = 0;
		objPoint(1) = -line(2) / line(1);
	}
}


////////////////////////  the mathods to determin dmgage value

void
ovalConcrete::getDltDmgTensile(const Vector line0, const Vector lineu, Vector sig,
	double& dDmgt, double& dDmgs)
{

	Vector direction(3), origin(2);

	direction = getLineByTwoPoints(origin, sig);
	Vector crossP0(2), crossPu(2);

	getTwoLinesCrossPoint(crossP0, direction, line0);
	getTwoLinesCrossPoint(crossPu, direction, lineu);

	dDmgt = (sig - crossP0)(0) / (crossPu - crossP0).Norm();
	dDmgs = abs((sig - crossP0)(1)) / (crossPu - crossP0).Norm();

}

double
ovalConcrete::getDltDmgShear(const Vector line0, const Vector lineu, Vector sig)
{
	Vector direction(3);
	direction = getVerticalLine(sig);

	Vector crossP0(2), crossPu(2);

	getTwoLinesCrossPoint(crossP0, direction, line0);
	getTwoLinesCrossPoint(crossPu, direction, lineu);

	double dltDmg = (sig - crossP0).Norm() / (crossPu - crossP0).Norm();

	return dltDmg;
}

double
ovalConcrete::getDltDmgCompressive(const Vector ellipse0P, const Vector ellipseiP, Vector sig)
{
	Vector p(2), direction(3);
	p(0) = ellipse0P(2);
	direction = getLineByTwoPoints(sig, p);

	Vector crossP01(2), crossP02(2), crossPu1(2), crossPu2(2);

	getEllipseLineCrossPoint(crossP01, crossP02, ellipse0P, direction);
	getEllipseLineCrossPoint(crossPu1, crossPu2, ellipseiP, direction);

	Vector crossP0(2), crossPu(2);

	if (((sig - p) ^ (crossP01 - p)) > 0)
		crossP0 = crossP01;
	else
		crossP0 = crossP02;

	if (((sig - p) ^ (crossPu1 - p)) > 0)
		crossPu = crossPu1;
	else
		crossPu = crossPu2;

	double dltDmgc = (sig - crossP0).Norm() / (crossPu - crossP0).Norm();

	return dltDmgc;
}


////////////////////////  the mathods to calculate stress

void
ovalConcrete::getStressTensile(Vector sig, Vector& sigCor, double dmgt, double dmgs)
{
	Vector origin(2), direction(3);
	direction = getLineByTwoPoints(sig, origin);

	Vector p1(2), p2(2), line(3);
	if (sig(1) > 0)
	{
		p1(0) = (1 - dmgt) * ft0init;
		p2(1) = (1 - dmgs) * fs0init;
		line = getLineByTwoPoints(p1, p2);
	}
	else
	{
		p1(0) = (1 - dmgt) * ft0init;
		p2(1) = -(1 - dmgs) * fs0init;
		line = getLineByTwoPoints(p1, p2);
	}

	if (dmgt < 1.0 - 1e-10)
		getTwoLinesCrossPoint(sigCor, direction, line);
	else
		sigCor.Zero();
}

void
ovalConcrete::getStressShear(Vector sig, Vector& sigCor, double dmgc, double dmgts, double dmgcs)
{
	Vector direction(3);
	direction = getVerticalLine(sig);

	Vector p1(2), p2(2), line(3);

	double axRat = axRatMax * (1 - dmgcs) + axRatMin * dmgcs;
	if (sig(1) > 0)
	{
		p1(0) = (1 - dmgc) * fc0init / 2;
		p1(1) = abs(p1(0)) * axRat;
		p2(1) = (1 - dmgts) * fs0init;
		line = getLineByTwoPoints(p1, p2);
	}
	else
	{
		p1(0) = (1 - dmgc) * fc0init / 2;
		p1(1) = -abs(p1(0)) * axRat;
		p2(1) = -(1 - dmgts) * fs0init;
		line = getLineByTwoPoints(p1, p2);
	}
	getTwoLinesCrossPoint(sigCor, direction, line);
}


void
ovalConcrete::getStressCompressive(Vector sig, Vector& sigCor, double dmgc, double dmgcsc)
{
	Vector p(2), direction(3);
	p(0) = (1 - dmgc) * fc0init / 2;
	direction = getLineByTwoPoints(sig, p);

	double axRat = axRatMax * (1 - dmgcsc) + axRatMin * dmgcsc;

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
ovalConcrete::getTangentTensile(Vector sig, Vector sigCor, const Vector line0P, const Vector lineiP)
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

	double dDmgtPtb, dDmgsPtb, dmgtPtb, dmgsPtb;
	Vector sigTrCor(2);

	getDltDmgTensile(line0P, lineiP, sigTrN, dDmgtPtb, dDmgsPtb);
	dmgtPtb = dmgtP + (1 - dmgtP) * dDmgtPtb;
	if (sigTrN(1) > 0)
		dmgsPtb = dmgtspP + (1 - dmgtspP) * dDmgsPtb;
	else
		dmgsPtb = dmgtsnP + (1 - dmgtsnP) * dDmgsPtb;
	getStressTensile(sigTrN, sigTrCor, dmgtPtb, dmgsPtb);
	if (dmgtPtb > 1.0)
		dmgtPtb = 1.0;
	if (dmgsPtb > 1.0)
		dmgsPtb = 1.0;
	tan(0, 0) = (sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = (sigTrCor(1) - sigCor(1)) / perturb;

	getDltDmgTensile(line0P, lineiP, sigTrS, dDmgtPtb, dDmgsPtb);
	dmgtPtb = dmgtP + (1 - dmgtP) * dDmgtPtb;
	if (sigTrS(1) > 0)
		dmgsPtb = dmgtspP + (1 - dmgtspP) * dDmgsPtb;
	else
		dmgsPtb = dmgtsnP + (1 - dmgtsnP) * dDmgsPtb;
	if (dmgtPtb > 1.0)
		dmgtPtb = 1.0;
	if (dmgsPtb > 1.0)
		dmgsPtb = 1.0;
	getStressTensile(sigTrS, sigTrCor, dmgtPtb, dmgsPtb);
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


Matrix
ovalConcrete::getTangentShear(Vector sig, Vector sigCor, const Vector line0P, const Vector lineiP)
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

	double dDmgsPtb, dDmgtPtb, dmgtsPtb, dmgcsPtb, dmgtPtb;
	Vector sigTrCor(2);

	if (sigTrN(0) > 0)
	{
		Vector l0(3), li(3), p1(2), p2(2);
		p1(0) = ft0P;
		if (sigTrN(1) > 0)
		{
			p2(1) = fsp0P;
			l0 = getLineByTwoPoints(p1, p2);
		}
		else
		{
			p2(1) = -fsn0P;
			l0 = getLineByTwoPoints(p1, p2);
		}
		p1(0) = ftiP;
		if (sigTrN(1) > 0)
		{
			p2(1) = ftspiP;
			li = getLineByTwoPoints(p1, p2);
		}
		else
		{
			p2(1) = -ftsniP;
			li = getLineByTwoPoints(p1, p2);
		}
		getDltDmgTensile(l0, li, sigTrN, dDmgtPtb, dDmgsPtb);
		dmgtPtb = dmgtP + (1 - dmgtP) * dDmgtPtb;
		dmgtsPtb = dmgtspP + (1 - dmgtspP) * dDmgsPtb;
		dmgcsPtb = dmgcspP;
		getStressTensile(sigTrN, sigTrCor, dmgtPtb, dmgtsPtb);
	}
	else
	{
		dmgtPtb = dmgtP;
		dDmgsPtb = getDltDmgShear(line0P, lineiP, sigTrN);
		if (sig(1) > 0)
		{
			dmgtsPtb = dmgtspP + (1 - dmgtspP) * dDmgsPtb * sigTrN(0) / fc0P / 2;
			dmgcsPtb = dmgcspP + (1 - dmgcspP) * dDmgsPtb;
		}
		else
		{
			dmgtsPtb = dmgtsnP + (1 - dmgtsnP) * dDmgsPtb * sigTrN(0) / fc0P / 2;
			dmgcsPtb = dmgcsnP + (1 - dmgcsnP) * dDmgsPtb;
		}
		getStressShear(sigTrN, sigTrCor, dmgcP, dmgtsPtb, dmgcsPtb);
	}
	if (dmgtsPtb > 1.0)
		dmgtsPtb = 1.0;
	if (dmgcsPtb > ultmDmg)
		dmgcsPtb = ultmDmg;
	tan(0, 0) = (sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = (sigTrCor(1) - sigCor(1)) / perturb;


	dDmgsPtb = getDltDmgShear(line0P, lineiP, sigTrS);
	if (sig(1) > 0)
	{
		dmgtsPtb = dmgtspP + (1 - dmgtspP) * dDmgsPtb * sigTrS(0) / fc0P / 2;
		dmgcsPtb = dmgcspP + dDmgsPtb;
	}
	else
	{
		dmgtsPtb = dmgtsnP + (1 - dmgtsnP) * dDmgsPtb * sigTrS(0) / fc0P / 2;
		dmgcsPtb = dmgcsnP + dDmgsPtb;
	}
	if (dmgtsPtb > 1.0)
		dmgtsPtb = 1.0;
	if (dmgcsPtb > ultmDmg)
		dmgcsPtb = ultmDmg;
	getStressShear(sigTrS, sigTrCor, dmgcP, dmgtsPtb, dmgcsPtb);
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
ovalConcrete::getTangentShear(Vector sig, Vector sigCor, const Vector line0P)
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
ovalConcrete::getTangentCopressive(Vector sig, Vector sigCor, const Vector ellipse0P, const Vector ellipseiP)
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

	double dDmgp, dmgcp, dmgcsp;
	Vector sigTrCor(2);

	dDmgp = getDltDmgCompressive(ellipse0P, ellipseiP, sigTrN);
	dmgcp = dmgcP + (1 - dmgcP) * dDmgp;
	if (sig(1) > 0)
		dmgcsp = dmgcspP + (1 - dmgcspP) * dDmgp;
	else
		dmgcsp = dmgcsnP + (1 - dmgcsnP) * dDmgp;
	if (dmgcp > ultmDmg)
		dmgcp = ultmDmg;
	if (dmgcsp > ultmDmg)
		dmgcsp = ultmDmg;
	getStressCompressive(sigTrN, sigTrCor, dmgcp, dmgcsp);
	tan(0, 0) = -(sigTrCor(0) - sigCor(0)) / perturb;
	tan(1, 0) = -(sigTrCor(1) - sigCor(1)) / perturb;

	dDmgp = getDltDmgCompressive(ellipse0P, ellipseiP, sigTrS);
	dmgcp = dmgcP + (1 - dmgcP) * dDmgp;
	if (sig(1) > 0)
		dmgcsp = dmgcspP + (1 - dmgcspP) * dDmgp;
	else
		dmgcsp = dmgcsnP + (1 - dmgcsnP) * dDmgp;
	if (dmgcp > ultmDmg)
		dmgcp = ultmDmg;
	if (dmgcsp > ultmDmg)
		dmgcsp = ultmDmg;
	getStressCompressive(sigTrS, sigTrCor, dmgcp, dmgcsp);
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
ovalConcrete::getTangentCopressive(Vector sig, Vector sigCor, const Vector ellipse0P)
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
ovalConcrete::getTrialStressTangent(Vector& st, Matrix& tgt, Vector et)
{
	st = sigP + eP * (et - epsP);
	tgt = eP;

	if (et(0) >= epsresn)
	{
		double sigMax, sigMin;
		sigMax = sigP(0) + knInit * (et(0) - epsP(0));
		sigMin = ktP * (et(0) - epsresn);

		if (st(0) > sigMin)
		{
			if (st(0) > sigMax)
			{
				st(0) = sigMax;
				tgt(0, 0) = knInit;
				if (st(0) < sigMin)
				{
					st(0) = sigMin;
					tgt(0, 0) = ktP;
				}
			}
		}

		if (eps(1) >= epsressp)
		{
			double tauMax, tauMin;
			tauMax = sigP(1) + ksInit * (et(1) - epsP(1));
			tauMin = ktspP * (et(1) - epsressp);
			if (st(1) > tauMin)
			{
				if (st(1) > tauMax)
				{
					st(1) = tauMax;
					tgt(1, 1) = ksInit;
					if (st(1) < tauMin)
					{
						st(1) = tauMin;
						tgt(1, 1) = ktspP;
					}
				}
			}
		}
		else if (eps(1) <= epsressn)
		{
			double tauMax, tauMin;
			tauMax = sigP(1) + ksInit * (et(1) - epsP(1));
			tauMin = ktsnP * (et(1) - epsressn);
			if (st(1) < tauMin)
			{
				if (st(1) < tauMax)
				{
					st(1) = tauMax;
					tgt(1, 1) = ksInit;
					if (st(1) > tauMin)
					{
						st(1) = tauMin;
						tgt(1, 1) = ktsnP;
					}
				}
			}
		}
		else
		{
			st(1) = 0.0;
			tgt(1, 1) = 1.0e-10;
		}
	}
	else if (et(0) <= epsresn)
	{
		double sigMax, sigMin;
		sigMax = sigP(0) + knInit * (et(0) - epsP(0));
		sigMin = kcP * (et(0) - epsresn);
		if (st(0) < sigMin)
		{
			if (st(0) < sigMax)
			{
				st(0) = sigMax;
				tgt(0, 0) = knInit;
				if (st(0) > sigMin)
				{
					st(0) = sigMin;
					tgt(0, 0) = kcP;
				}
			}
		}

		if (eps(1) >= epsressp)
		{
			double tauMax, tauMin;
			tauMax = sigP(1) + ksInit * (et(1) - epsP(1));
			tauMin = kcspP * (et(1) - epsressp);
			if (st(1) > tauMin)
			{
				if (st(1) > tauMax)
				{
					st(1) = tauMax;
					tgt(1, 1) = ksInit;
					if (st(1) < tauMin)
					{
						st(1) = tauMin;
						tgt(1, 1) = kcspP;
					}
				}
			}
		}
		else if (eps(1) <= epsressn)
		{
			double tauMax, tauMin;
			tauMax = sigP(1) + ksInit * (et(1) - epsP(1));
			tauMin = kcsnP * (et(1) - epsressn);
			if (st(1) < tauMin)
			{
				if (st(1) < tauMax)
				{
					st(1) = tauMax;
					tgt(1, 1) = ksInit;
					if (st(1) > tauMin)
					{
						st(1) = tauMin;
						tgt(1, 1) = kcsnP;
					}
				}
			}
		}
		else
		{
			st(1) = 0.0;
			tgt(1, 1) = 1.0e-10;
		}
	}
}