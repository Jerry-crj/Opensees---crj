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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionTimoQualShear.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSectionTimoQualShear.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionTimoQualShear.h>
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

ID FiberSectionTimoQualShear::code(8);

void 
FiberSectionTimoQualShear::setSectionSize(double h0, double b0)
{
	b = b0;
	h = h0;
}


int
FiberSectionTimoQualShear::addFiber(Fiber &newFiber)
{
	// need to create a larger array
	if (numFibers == sizeFibers) {
		int newSize = 2 * sizeFibers;
		NDMaterial **newArray = new NDMaterial *[newSize];
		double *newMatData = new double[3 * newSize];

		if (newArray == 0 || newMatData == 0) {
			opserr << "FiberSectionTimoQualShear::addFiber -- failed to allocate Fiber pointers\n";
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
		opserr << "FiberSectionTimoQualShear::addFiber -- failed to get copy of a Material\n";
		return -1;
	}

	numFibers++;

	// Recompute centroid
	Abar += Area;
	QzBar += yLoc * Area;

	yBar = QzBar / Abar;

	return 0;
}

FiberSectionTimoQualShear::FiberSectionTimoQualShear(int tag, int num, Fiber** fibers) :
	SectionForceDeformation(tag, SEC_TAG_FiberSectionTimoQualShear),
	numFibers(num), sizeFibers(num), theMaterials(0), matData(0), b(0), h(0), 
	QzBar(0.0), Abar(0.0), yBar(0.0), sectionIntegr(0), e(4), s(4), ks(4, 4), ki(0)
{
	if (numFibers != 0) {
		theMaterials = new NDMaterial *[numFibers];

		if (theMaterials == 0) {
			opserr << "FiberSectionTimoQualShearN::FiberSectionTimoQualShearN -- failed to allocate Material pointers\n";
			exit(-1);
		}

		matData = new double[numFibers * 3];

		if (matData == 0) {
			opserr << "FiberSectionTimoQualShearN::FiberSectionTimoQualShearN -- failed to allocate double array for material data\n";
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
				opserr << "FiberSectionTimoQualShearN::FiberSectionTimoQualShearN -- failed to get copy of a Material\n";
				exit(-1);
			}
		}

		yBar = QzBar / Abar;
	}

}


// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionTimoQualShear::FiberSectionTimoQualShear() :
	SectionForceDeformation(0, SEC_TAG_FiberSectionTimoQualShear),
	numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
	QzBar(0.0), Abar(0.0), yBar(0.0), b(0), h(0),
	sectionIntegr(0), e(4), s(4), ks(4, 4), ki(0)
{

}

// destructor:
FiberSectionTimoQualShear::~FiberSectionTimoQualShear()
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

	if (ki != 0)
		delete ki;
}

int
FiberSectionTimoQualShear::setTrialSectionDeformation(const Vector &deforms)
{
	if (b == 0 || h == 0)
		opserr << "the section size has not been specified " << endln;

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

		for (int i = 0; i < numFibers; i++)
		{
			yLocs[i] = matData[3 * i];
			zLocs[i] = matData[3 * i + 1];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 4), crdCoeffTran(4, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = -y;
		crdCoeff(0, 2) = -4. / 3. * pow(y, 3) / pow(h, 2);
		crdCoeff(1, 3) = 1 - pow(2 * y / h, 2);

		crdCoeffTran.addMatrixTranspose(0.0, crdCoeff, 1.0);

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
FiberSectionTimoQualShear::getInitialTangent(void)
{
	if (ki != 0)
		return *ki;
	
	ki = new Matrix(4, 4);

	static double yLocs[10000];
	static double zLocs[10000];
	static double fiberArea[10000];

	if (sectionIntegr != 0) 
	{
		sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
		sectionIntegr->getFiberWeights(numFibers, fiberArea);
	}
	else
	{
		for (int i = 0; i < numFibers; i++) 
		{
			yLocs[i] = matData[3 * i];
			zLocs[i] = matData[3 * i + 1];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) 
	{
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double A = fiberArea[i];

		Matrix crdCoeff(2, 4), crdCoeffTran(4, 2);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = -y;
		crdCoeff(0, 2) = -4 * pow(y, 3) / 3 / pow(h, 2);
		crdCoeff(1, 3) = 1 - pow(2 * y / h, 2);

		crdCoeffTran.addMatrixTranspose(0.0, crdCoeff, 1.0);

		const Matrix &fiberTangent = theMat->getInitialTangent();

		*ki += (crdCoeffTran * fiberTangent * crdCoeff) * A;
	}

	return *ki;
}

const Vector&
FiberSectionTimoQualShear::getSectionDeformation(void)
{
	return e;
}

const Matrix&
FiberSectionTimoQualShear::getSectionTangent(void)
{
	return ks;
}

const Vector&
FiberSectionTimoQualShear::getStressResultant(void)
{
	return s;
}

SectionForceDeformation*
FiberSectionTimoQualShear::getCopy(void)
{
	FiberSectionTimoQualShear *theCopy = new FiberSectionTimoQualShear();
	theCopy->setTag(this->getTag());

	theCopy->numFibers = numFibers;
	theCopy->sizeFibers = numFibers;
	theCopy->h = h;
	theCopy->b = b;

	if (numFibers != 0) {
		theCopy->theMaterials = new NDMaterial *[numFibers];

		if (theCopy->theMaterials == 0) {
			opserr << "FiberSectionTimoQualShear::FiberSectionTimoQualShear -- failed to allocate Material pointers\n";
			exit(-1);
		}

		theCopy->matData = new double[numFibers * 3];

		if (theCopy->matData == 0) {
			opserr << "FiberSectionTimoQualShear::FiberSectionTimoQualShear -- failed to allocate double array for material data\n";
			exit(-1);
		}


		for (int i = 0; i < numFibers; i++)
		{
			theCopy->matData[i * 3] = matData[i * 3];
			theCopy->matData[i * 3 + 1] = matData[i * 3 + 1];
			theCopy->matData[i * 3 + 2] = matData[i * 3 + 2];
			theCopy->theMaterials[i] = theMaterials[i]->getCopy();

			if (theCopy->theMaterials[i] == 0) {
				opserr << "FiberSectionTimoQualShear::getCopy -- failed to get copy of a Material\n";
				exit(-1);
			}
		}
	}

	theCopy->e = e;
	theCopy->QzBar = QzBar;
	theCopy->Abar = Abar;
	theCopy->yBar = yBar;

	theCopy->ks = ks;
	theCopy->s = s;

	if (sectionIntegr != 0)
		theCopy->sectionIntegr = sectionIntegr->getCopy();

	if (ki != 0)
		theCopy->ki = ki;

	return theCopy;
}

const ID&
FiberSectionTimoQualShear::getType()
{
	return code;
}

int
FiberSectionTimoQualShear::getOrder() const
{
	return 8;
}

int
FiberSectionTimoQualShear::commitState(void)
{	
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->commitState();

	return err;
}

int
FiberSectionTimoQualShear::revertToLastCommit(void)
{
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->revertToLastCommit();

	return err;
}

int
FiberSectionTimoQualShear::revertToStart(void)
{
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->revertToStart();

	return err;
}

int
FiberSectionTimoQualShear::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	return res;
}

int
FiberSectionTimoQualShear::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	return res;
}

void
FiberSectionTimoQualShear::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "\tSection type: FiberSectionTimoQualShear, tag: " << this->getTag() << endln;
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
FiberSectionTimoQualShear::setResponse(const char **argv, int argc, OPS_Stream &output)
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
	else if (strcmp(argv[0], "damgaeFactor") == 0 || strcmp(argv[0], "dam") == 0)
	{
		Vector res(2 * numFibers);
		theResponse = new MaterialResponse(this, 3, res);
	}
	else if (strcmp(argv[0], "secForce") == 0)
	{
		Vector res(3);
		theResponse = new MaterialResponse(this, 4, res);
	}
	else if (strcmp(argv[0], "secDeform") == 0)
	{
		Vector res(3);
		theResponse = new MaterialResponse(this, 5, res);
	}
	else
		opserr << "need to be add" << endln;

	return theResponse;
}


int
FiberSectionTimoQualShear::getResponse(int responseID, Information &sectInfo)
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
		Vector data(numFibers * 2);
		int count = 0;
		for (int i = 0; i < numFibers; i++)
		{
			ASIConcrete* mat;
			if (strcmp(theMaterials[i]->getType(), "ASIConcrete") == 0)
			{
				mat = (ASIConcrete*)theMaterials[i];

				Vector res(2);
				//res = mat->getDamage();
				data[count] = res(0);
				data[count + 1] = res(1);
				count += 2;
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
		return sectInfo.setVector(e);
	}
	default:
		return -1;
	}
}