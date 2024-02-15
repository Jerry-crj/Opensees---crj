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

// $Revision: 1.32 $
// $Date: 2010-08-16 05:05:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionTimoConstShear.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSectionTimoConstShear.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionTimoConstShear.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>
#include <ElasticMaterial.h>
#include <SectionIntegration.h>
#include <elementAPI.h>
#include <iostream>
#include <string>
#include <direct.h>
#include <ASIConcrete.h>

ID FiberSectionTimoConstShear::code(4);

int
FiberSectionTimoConstShear::addFiber(Fiber &newFiber)
{
	// need to create a larger array
	if (numFibers == sizeFibers) {
		int newSize = 2 * sizeFibers;
		NDMaterial **newArray = new NDMaterial *[newSize];
		double *newMatData = new double[3 * newSize];

		if (newArray == 0 || newMatData == 0) {
			opserr << "FiberSectionTimoConstShear::addFiber -- failed to allocate Fiber pointers\n";
			exit(-1);
		}

		// copy the old pointers
		for (int i = 0; i < numFibers; i++) {
			newArray[i] = theMaterials[i];
			newMatData[3 * i] = matData[3 * i];
			newMatData[3 * i + 1] = matData[3 * i + 1];
			newMatData[3 * i + 2] = matData[3 * i + 2];
		}

		// initialize new memomry
		for (int i = numFibers; i < newSize; i++) {
			newArray[i] = 0;
			newMatData[3 * i] = 0.0;
			newMatData[3 * i + 1] = 0.0;
			newMatData[3 * i + 2] = 0.0;
		}
		sizeFibers = newSize;

		// set new memory
		if (theMaterials != 0) {
			delete[] theMaterials;
			delete[] matData;
		}

		theMaterials = newArray;
		matData = newMatData;
	}

	// set the new pointers
	double yLoc, zLoc, Area;
	newFiber.getFiberLocation(yLoc, zLoc);
	Area = newFiber.getArea();
	matData[numFibers * 3] = yLoc;
	matData[numFibers * 3 + 1] = zLoc;
	matData[numFibers * 3 + 2] = Area;
	NDMaterial *theMat = newFiber.getNDMaterial();
	theMaterials[numFibers] = theMat->getCopy();

	if (theMaterials[numFibers] == 0) {
		opserr << "FiberSectionTimoConstShear::addFiber -- failed to get copy of a Material\n";
		return -1;
	}

	numFibers++;

	// Recompute centroid
	Abar += Area;
	QzBar += yLoc * Area;

	yBar = QzBar / Abar;

	return 0;
}

FiberSectionTimoConstShear::FiberSectionTimoConstShear(int tag, int num, Fiber **fibers) :
	SectionForceDeformation(tag, SEC_TAG_FiberSectionTimoConstShear),
	numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
	QzBar(0.0), Abar(0.0), yBar(0.0), sectionIntegr(0), e(3), s(3), ks(3, 3), reduce(2, 2), kInitial(3, 3), eP(3)
{
	if (numFibers != 0) {
		theMaterials = new NDMaterial *[numFibers];

		if (theMaterials == 0) {
			opserr << "FiberSectionTimoConstShearN::FiberSectionTimoConstShearN -- failed to allocate Material pointers\n";
			exit(-1);
		}

		matData = new double[numFibers * 3];

		if (matData == 0) {
			opserr << "FiberSectionTimoConstShearN::FiberSectionTimoConstShearN -- failed to allocate double array for material data\n";
			exit(-1);
		}

		for (int i = 0; i < numFibers; i++) {
			Fiber *theFiber = fibers[i];
			double yLoc, zLoc, Area;
			theFiber->getFiberLocation(yLoc, zLoc);
			Area = theFiber->getArea();

			QzBar += yLoc * Area;
			Abar += Area;

			matData[i * 3] = yLoc;
			matData[i * 3 + 1] = zLoc;
			matData[i * 3 + 2] = Area;
			NDMaterial *theMat = theFiber->getNDMaterial();
			theMaterials[i] = theMat->getCopy();

			if (theMaterials[i] == 0) {
				opserr << "FiberSectionTimoConstShearN::FiberSectionTimoConstShearN -- failed to get copy of a Material\n";
				exit(-1);
			}
		}

		yBar = QzBar / Abar;
	}

	//s = new Vector(sData, 3);
	//ks = new Matrix(kData, 3, 3);

	for (int i = 0; i < 3; i++)
		sData[i] = 0.0;

	for (int i = 0; i < 9; i++)
		kData[i] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_MY;
	code(3) = SECTION_RESPONSE_T;

}


// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionTimoConstShear::FiberSectionTimoConstShear() :
	SectionForceDeformation(0, SEC_TAG_FiberSectionTimoConstShear),
	numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
	QzBar(0.0), Abar(0.0), yBar(0.0), sectionIntegr(0), e(3), s(3), ks(3,3)
{
	//s = new Vector(sData, 3);
	//ks = new Matrix(kData, 3, 3);

	for (int i = 0; i < 3; i++)
		sData[i] = 0.0;

	for (int i = 0; i < 9; i++)
		kData[i] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_MY;
	code(3) = SECTION_RESPONSE_T;
}

// destructor:
FiberSectionTimoConstShear::~FiberSectionTimoConstShear()
{
	if (theMaterials != 0) {
		for (int i = 0; i < numFibers; i++)
			if (theMaterials[i] != 0)
				delete theMaterials[i];

		delete[] theMaterials;
	}

	if (matData != 0)
		delete[] matData;

	if (sectionIntegr != 0)
		delete sectionIntegr;

}



int
FiberSectionTimoConstShear::setTrialSectionDeformation(const Vector &deforms)
{
	int err = 0;
	e = deforms;

	s.Zero();
	ks.Zero();

	static double yLocs[10000];
	static double zLocs[10000];
	static double fiberArea[10000];

	if (sectionIntegr != 0) {
		sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
		sectionIntegr->getFiberWeights(numFibers, fiberArea);
	}
	else {

		for (int i = 0; i < numFibers; i++) {

			yLocs[i] = matData[3 * i];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 3), crdCoeffTran(3, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(1, 2) = 1;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
				crdCoeffTran(j, i) = crdCoeff(i, j);

		Vector fiberStrain(2);
		fiberStrain.addMatrixVector(1.0, crdCoeff, e, 1.0);
		err += theMat->setTrialStrain(fiberStrain);

		const Vector &fiberStress = theMat->getStress();
		const Matrix &fiberTangent = theMat->getTangent();

		s += (crdCoeffTran * fiberStress) * A;
		ks += (crdCoeffTran * fiberTangent * crdCoeff) * A;
	}

	return err;
}

const Matrix&
FiberSectionTimoConstShear::getInitialTangent(void)
{
	kInitial.Zero();

	static double yLocs[10000];
	static double zLocs[10000];
	static double fiberArea[10000];

	if (sectionIntegr != 0) {
		sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
		sectionIntegr->getFiberWeights(numFibers, fiberArea);
	}
	else {

		for (int i = 0; i < numFibers; i++) {

			yLocs[i] = matData[3 * i];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 3), crdCoeffTran(3, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(1, 2) = 1;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
				crdCoeffTran(j, i) = crdCoeff(i, j);

		const Matrix &fiberTangent = theMat->getInitialTangent();

		kInitial += (crdCoeffTran * fiberTangent * crdCoeff) * A;
	}

	return kInitial;
}

const Vector&
FiberSectionTimoConstShear::getSectionDeformation(void)
{
	return e;
}

const Matrix&
FiberSectionTimoConstShear::getSectionTangent(void)
{
	return ks;
}

const Vector&
FiberSectionTimoConstShear::getStressResultant(void)
{
	return s;
}

SectionForceDeformation*
FiberSectionTimoConstShear::getCopy(void)
{
	FiberSectionTimoConstShear *theCopy = new FiberSectionTimoConstShear();
	theCopy->setTag(this->getTag());

	theCopy->numFibers = numFibers;
	theCopy->sizeFibers = numFibers;

	if (numFibers != 0) {
		theCopy->theMaterials = new NDMaterial *[numFibers];

		if (theCopy->theMaterials == 0) {
			opserr << "FiberSectionTimoConstShear::FiberSectionTimoConstShear -- failed to allocate Material pointers\n";
			exit(-1);
		}

		theCopy->matData = new double[numFibers * 3];

		if (theCopy->matData == 0) {
			opserr << "FiberSectionTimoConstShear::FiberSectionTimoConstShear -- failed to allocate double array for material data\n";
			exit(-1);
		}


		for (int i = 0; i < numFibers; i++) {
			theCopy->matData[i * 3] = matData[i * 3];
			theCopy->matData[i * 3 + 1] = matData[i * 3 + 1];
			theCopy->matData[i * 3 + 2] = matData[i * 3 + 2];
			theCopy->theMaterials[i] = theMaterials[i]->getCopy();

			if (theCopy->theMaterials[i] == 0) {
				opserr << "FiberSectionTimoConstShear::getCopy -- failed to get copy of a Material\n";
				exit(-1);
			}
		}
	}

	theCopy->e = e;
	theCopy->QzBar = QzBar;
	theCopy->Abar = Abar;
	theCopy->yBar = yBar;

	for (int i = 0; i < 9; i++)
		theCopy->kData[i] = kData[i];

	for (int i = 0; i < 3; i++)
		theCopy->sData[i] = kData[i];

	if (sectionIntegr != 0)
		theCopy->sectionIntegr = sectionIntegr->getCopy();
	else
		theCopy->sectionIntegr = 0;

	theCopy->reduce = reduce;
	theCopy->kInitial = kInitial;
	theCopy->eP = eP;

	return theCopy;
}

const ID&
FiberSectionTimoConstShear::getType()
{
	return code;
}

int
FiberSectionTimoConstShear::getOrder() const
{
	return 4;
}

int
FiberSectionTimoConstShear::commitState(void)
{	
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->commitState();

	return err;
}

int
FiberSectionTimoConstShear::revertToLastCommit(void)
{
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->revertToLastCommit();

	return err;
}

int
FiberSectionTimoConstShear::revertToStart(void)
{
	// revert the fibers to start    
	int err = 0;

	s.Zero();
	ks.Zero();

	static double yLocs[10000];
	static double zLocs[10000];
	static double fiberArea[10000];

	if (sectionIntegr != 0) {
		sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
		sectionIntegr->getFiberWeights(numFibers, fiberArea);
	}
	else {
		for (int i = 0; i < numFibers; i++) {
			yLocs[i] = matData[3 * i];
			fiberArea[i] = matData[3 * i + 1];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		// invoke revertToStart on the material
		// invoke revertToLast on the material
		err += theMat->revertToStart();

		const Vector &stress = theMat->getStress();
		const Matrix &tangent = theMat->getTangent();

		Matrix crdCoeff(2, 3);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(1, 2) = 1;

		Matrix secStiff(3, 3), CoeffD(2, 3);
		Vector secForce(3);

		secForce.addMatrixTransposeVector(1.0, crdCoeff, stress, 1.0);
		CoeffD.addMatrixProduct(0.0, tangent, crdCoeff, 1.0);
		secStiff.addMatrixTransposeProduct(0.0, crdCoeff, CoeffD, 1.0);
		s += secForce * A;
		ks += secStiff * A;
	}

	return err;
}

int
FiberSectionTimoConstShear::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	return res;
}

int
FiberSectionTimoConstShear::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	return res;
}

void
FiberSectionTimoConstShear::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "\tSection type: FiberSectionTimoConstShear, tag: " << this->getTag() << endln;
		s << "\tNumber of Fibers: " << numFibers << endln;
		s << "\tCentroid: (" << -yBar << ", " << ')' << endln;

		s << "\nMaterial Information :";
		if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
			for (int i = 0; i < numFibers; i++) {
				s << "\n\tFiber  Tag: " << i + 1;
				s << "\n\tLocation y = " << matData[2 * i] << "\tArea = " << matData[2 * i + 1] << endln;
				theMaterials[i]->Print(s, flag);
			}
		}
	}
}

Response*
FiberSectionTimoConstShear::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	if (argc > 2 || strcmp(argv[0], "fiber") == 0) {

		int key = numFibers;
		int passarg = 2;

		if (argc == 4)   // fiber number was input directly
		{
			double yCoord = atof(argv[1]);
			double zCoord = atof(argv[2]);
			double closestDist;
			double ySearch, zSearch, dy, dz;
			double distance;
			ySearch = matData[0];
			zSearch = matData[1];
			dy = ySearch - yCoord;
			dz = zSearch - zCoord;
			closestDist = sqrt(dy * dy + dz * dz);
			key = 0;
			for (int j = 1; j < numFibers; j++) {
				ySearch = matData[3 * j];
				zSearch = matData[3 * j + 1];
				dy = ySearch - yCoord;
				dz = zSearch - zCoord;
				distance = sqrt(dy * dy + dz * dz);
				if (distance < closestDist) {
					closestDist = distance;
					key = j;
				}
			}
			passarg = 3;
			if (key < numFibers && key >= 0)
				theResponse = theMaterials[key]->setResponse(&argv[passarg], argc - passarg, output);
		}
	}
	else if (strcmp(argv[0], "strain") == 0)
	{
		Vector res(2 * numFibers);
		theResponse = new MaterialResponse(this, 1, res);
	}
	else if (strcmp(argv[0], "stress") == 0)
	{
		Vector res(2 * numFibers);
		theResponse = new MaterialResponse(this, 2, res);
	}
	else if (strcmp(argv[0], "damgae") == 0 || strcmp(argv[0], "dmg") == 0)
	{
		Vector res(3 * numFibers);
		theResponse = new MaterialResponse(this, 3, res);
	}
	else if (strcmp(argv[0], "tangent") == 0)
	{
		Vector res(6 * numFibers);
		theResponse = new MaterialResponse(this, 4, res);
	}
	else if (strcmp(argv[0], "paperInfo") == 0)
	{
		Vector res(7 * numFibers);
		theResponse = new MaterialResponse(this, 5, res);
	}
	else
		opserr << "need to be added" << endln;

	return theResponse;
}


int
FiberSectionTimoConstShear::getResponse(int responseID, Information &sectInfo)
{
	switch (responseID) 
	{
	case 1:
	{
		Vector data(2 * numFibers);
		for (int i = 0; i < numFibers; i++)
		{
			Vector matStrain(2);
			matStrain = theMaterials[i]->getStrain();
			data(2 * i) = matStrain(0);
			data(2 * i + 1) = matStrain(1);
		}
		return sectInfo.setVector(data);
	}
	case 2:
	{
		Vector data(2 * numFibers);
		for (int i = 0; i < numFibers; i++)
		{
			Vector matStress(2);
			matStress = theMaterials[i]->getStress();
			data(2 * i) = matStress(0);
			data(2 * i + 1) = matStress(1);
		}
		return sectInfo.setVector(data);
	}
	case 3:
	{
		Vector data(numFibers * 3);
		int count = 0;
		for (int i = 0; i < numFibers; i++)
		{
			ASIConcrete* theMat;
			if (strcmp(theMaterials[i]->getType(), "ASIConcrete") == 0)
			{
				theMat = (ASIConcrete*)theMaterials[i];

				Vector res = theMat->getDamage();
				data[count + 0] = res(0);
				data[count + 1] = res(1);
				data[count + 2] = res(2);
				count += 3;
			}
		}
		return sectInfo.setVector(data);
	}
	case 4:
	{
		return sectInfo.setVector(s);
	}
	case 5:
	{
		Vector data(7 * numFibers);
		for (int i = 0; i < numFibers; i++)
		{
			data(7 * i + 0) = theMaterials[i]->getTag();
			data(7 * i + 1) = matData[i * 3];
			data(7 * i + 2) = matData[i * 3 + 1];
			data(7 * i + 3) = theMaterials[i]->getStrain()[0];
			data(7 * i + 4) = theMaterials[i]->getStress()[0];
			data(7 * i + 5) = theMaterials[i]->getStrain()[1];
			data(7 * i + 6) = theMaterials[i]->getStress()[1];
		}
		return sectInfo.setVector(data);
	}
	default:
		return -1;
	}
}


const Matrix
FiberSectionTimoConstShear::perturbToGetTangent(void)
{
	Matrix stiff(3, 3);

	double err = 0;

	static double yLocs[10000];
	static double zLocs[10000];
	static double fiberArea[10000];

	if (sectionIntegr != 0) {
		sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
		sectionIntegr->getFiberWeights(numFibers, fiberArea);
	}
	else {

		for (int i = 0; i < numFibers; i++) {

			yLocs[i] = matData[2 * i];
			fiberArea[i] = matData[2 * i + 1];
		}
	}

	double pertDisp = 1e-6;
	double pertAngle = 1e-8;

	Vector epsTr(3), sigTr, force, forceTr(3);

	force = s;

	forceTr.Zero();
	epsTr = e;
	if (e(0) - eP(0) >= 0)
		epsTr(0) += pertDisp;
	else
		epsTr(0) -= pertDisp;

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 3), crdCoeffTran(3, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(1, 2) = 1;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
				crdCoeffTran(j, i) = crdCoeff(i, j);

		Vector fiberStrainTr(2);
		fiberStrainTr.addMatrixVector(1.0, crdCoeff, epsTr, 1.0);
		err += theMat->setTrialStrain(fiberStrainTr);
		const Vector &fiberStressTr = theMat->getStress();

		forceTr += (crdCoeffTran * fiberStressTr) * A;
	}

	if (e(0) - eP(0) >= 0)
	{
		stiff(0, 0) = (forceTr(0) - force(0)) / pertDisp;
		stiff(1, 0) = (forceTr(1) - force(1)) / pertDisp;
		stiff(2, 0) = (forceTr(2) - force(2)) / pertDisp;
	}
	else
	{
		stiff(0, 0) = -(forceTr(0) - force(0)) / pertDisp;
		stiff(1, 0) = -(forceTr(1) - force(1)) / pertDisp;
		stiff(2, 0) = -(forceTr(2) - force(2)) / pertDisp;
	}

	forceTr.Zero();
	epsTr = e;
	if (e(1) - eP(1) >= 0)
		epsTr(1) += pertAngle;
	else
		epsTr(1) -= pertAngle;

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 3), crdCoeffTran(3, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(1, 2) = 1;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
				crdCoeffTran(j, i) = crdCoeff(i, j);

		Vector fiberStrainTr(2);
		fiberStrainTr.addMatrixVector(1.0, crdCoeff, epsTr, 1.0);
		err += theMat->setTrialStrain(fiberStrainTr);
		const Vector &fiberStressTr = theMat->getStress();

		forceTr += (crdCoeffTran * fiberStressTr) * A;
	}

	if (e(1) - eP(1) >= 0)
	{
		stiff(0, 1) = (forceTr(0) - force(0)) / pertAngle;
		stiff(1, 1) = (forceTr(1) - force(1)) / pertAngle;
		stiff(2, 1) = (forceTr(2) - force(2)) / pertAngle;
	}
	else
	{
		stiff(0, 1) = -(forceTr(0) - force(0)) / pertAngle;
		stiff(1, 1) = -(forceTr(1) - force(1)) / pertAngle;
		stiff(2, 1) = -(forceTr(2) - force(2)) / pertAngle;
	}

	forceTr.Zero();
	epsTr = e;
	if (e(2) - eP(2) >= 0)
		epsTr(2) += pertDisp;
	else
		epsTr(2) -= pertDisp;

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 3), crdCoeffTran(3, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(1, 2) = 1;

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
				crdCoeffTran(j, i) = crdCoeff(i, j);

		Vector fiberStrainTr(2);
		fiberStrainTr.addMatrixVector(1.0, crdCoeff, epsTr, 1.0);
		err += theMat->setTrialStrain(fiberStrainTr);
		const Vector &fiberStressTr = theMat->getStress();

		forceTr += (crdCoeffTran * fiberStressTr) * A;
	}

	if (e(2) - eP(2) >= 0)
	{
		stiff(0, 2) = (forceTr(0) - force(0)) / pertDisp;
		stiff(1, 2) = (forceTr(1) - force(1)) / pertDisp;
		stiff(2, 2) = (forceTr(2) - force(2)) / pertDisp;
	}
	else
	{
		stiff(0, 2) = -(forceTr(0) - force(0)) / pertDisp;
		stiff(1, 2) = -(forceTr(1) - force(1)) / pertDisp;
		stiff(2, 2) = -(forceTr(2) - force(2)) / pertDisp;
	}

	return stiff;
}