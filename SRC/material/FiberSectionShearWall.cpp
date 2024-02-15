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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionShearWall.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSectionShearWall.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionShearWall.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>
#include <ElasticMaterial.h>
#include <SectionIntegration.h>
#include <elementAPI.h>
#include <string.h>

ID FiberSectionShearWall::code(4);

int
FiberSectionShearWall::addFiber(Fiber &newFiber)
{
	// need to create a larger array
	if (numFibers == sizeFibers) {
		int newSize = 2 * sizeFibers;
		NDMaterial **newArray = new NDMaterial *[newSize];
		double *newMatData = new double[3 * newSize];

		if (newArray == 0 || newMatData == 0) {
			opserr << "FiberSectionShearWall::addFiber -- failed to allocate Fiber pointers\n";
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
		opserr << "FiberSectionShearWall::addFiber -- failed to get copy of a Material\n";
		return -1;
	}

	numFibers++;

	// Recompute centroid
	Abar += Area;
	QzBar += yLoc * Area;
	QyBar += zLoc * Area;

	yBar = QzBar / Abar;
	zBar = QyBar / Abar;

	return 0;
}

FiberSectionShearWall::FiberSectionShearWall(int tag, int num, Fiber **fibers) :
	SectionForceDeformation(tag, SEC_TAG_FiberSectionShearWall),
	numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
	QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), sectionIntegr(0), e(4), s(0), ks(0)
{
	if (numFibers != 0) {
		theMaterials = new NDMaterial *[numFibers];

		if (theMaterials == 0) {
			opserr << "FiberSectionShearWallN::FiberSectionShearWallN -- failed to allocate Material pointers\n";
			exit(-1);
		}

		matData = new double[numFibers * 3];

		if (matData == 0) {
			opserr << "FiberSectionShearWallN::FiberSectionShearWallN -- failed to allocate double array for material data\n";
			exit(-1);
		}

		for (int i = 0; i < numFibers; i++) {
			Fiber *theFiber = fibers[i];
			double yLoc, zLoc, Area;
			theFiber->getFiberLocation(yLoc, zLoc);
			Area = theFiber->getArea();

			QzBar += yLoc * Area;
			QyBar += zLoc * Area;
			Abar += Area;

			matData[i * 3] = yLoc;
			matData[i * 3 + 1] = zLoc;
			matData[i * 3 + 2] = Area;
			NDMaterial *theMat = theFiber->getNDMaterial();
			theMaterials[i] = theMat->getCopy();

			if (theMaterials[i] == 0) {
				opserr << "FiberSectionShearWallN::FiberSectionShearWallN -- failed to get copy of a Material\n";
				exit(-1);
			}
		}

		yBar = QzBar / Abar;
		zBar = QyBar / Abar;
	}

	s = new Vector(sData, 8);
	ks = new Matrix(kData, 8, 8);

	for (int i = 0; i < 8; i++)
		sData[i] = 0.0;

	for (int i = 0; i < 64; i++)
		kData[i] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_MY;
	code(3) = SECTION_RESPONSE_T;
}


// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionShearWall::FiberSectionShearWall() :
	SectionForceDeformation(0, SEC_TAG_FiberSectionShearWall),
	numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
	QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), sectionIntegr(0), e(4), s(0), ks(0)
{
	s = new Vector(sData, 8);
	ks = new Matrix(kData, 8, 8);

	for (int i = 0; i < 8; i++)
		sData[i] = 0.0;

	for (int i = 0; i < 64; i++)
		kData[i] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_MY;
	code(3) = SECTION_RESPONSE_T;
}

// destructor:
FiberSectionShearWall::~FiberSectionShearWall()
{
	if (theMaterials != 0) {
		for (int i = 0; i < numFibers; i++)
			if (theMaterials[i] != 0)
				delete theMaterials[i];

		delete[] theMaterials;
	}

	if (matData != 0)
		delete[] matData;

	if (s != 0)
		delete s;

	if (ks != 0)
		delete ks;

	if (sectionIntegr != 0)
		delete sectionIntegr;

}

int
FiberSectionShearWall::setTrialSectionDeformation(const Vector &deforms)
{
	int err = 0;
	e = deforms;

	s->Zero();
	ks->Zero();

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
			zLocs[i] = matData[3 * i + 1];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double z = zLocs[i] - zBar;
		double A = fiberArea[i];

		// determine material strain and set it
		Matrix crdCoeff(3, 8), crdCoeffTran(8, 3);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(0, 2) = -z;
		crdCoeff(0, 3) = -y * z;
		crdCoeff(1, 4) = 1;
		crdCoeff(1, 5) = -z;
		crdCoeff(2, 6) = 1;
		crdCoeff(2, 7) = -y;

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 8; j++)
				crdCoeffTran(j, i) = crdCoeff(i, j);

		Vector strain(3);
		strain.addMatrixVector(1.0, crdCoeff, e, 1.0);
		err += theMat->setTrialStrain(strain);

		const Vector &stress = theMat->getStress();
		const Matrix &tangent = theMat->getTangent();

		//if ( i < 60  && (abs(tangent(1, 1) - 1e6) > 1e-4 || abs(tangent(2, 2) - 1e6) > 1e-4) )
		//	double debug = 1;
		//else if ( i>=60 && (abs(tangent(1, 1)) > 1e-10 || abs(tangent(2, 2)) > 1e-10) )
		//	double debug = 1;

		Matrix secStiff(8, 8);
		Vector secForce(8);
		
		secForce = crdCoeffTran * stress;
		secStiff = crdCoeffTran * tangent * crdCoeff;

		(*s) += secForce * A;
		(*ks) += secStiff * A;
	}
	return err;
}

const Matrix&
FiberSectionShearWall::getInitialTangent(void)
{
	static double kInitialData[64];
	static Matrix kInitial(kInitialData, 8, 8);

	kInitial.Zero();

	static double yLocs[10000];
	static double zLocs[10000];
	static double fiberArea[10000];

	//if (sectionIntegr != 0) {
	//	sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
	//	sectionIntegr->getFiberWeights(numFibers, fiberArea);
	//}
	//else {
	//	for (int i = 0; i < numFibers; i++) {
	//		yLocs[i] = matData[3 * i];
	//		zLocs[i] = matData[3 * i + 1];
	//		fiberArea[i] = matData[3 * i + 2];
	//	}
	//}

	//for (int i = 0; i < numFibers; i++) {
	//	NDMaterial *theMat = theMaterials[i];
	//	double y = yLocs[i] - yBar;
	//	double z = zLocs[i] - zBar;
	//	double A = fiberArea[i];

	//	double tangent = theMat->getInitialTangent();

	//	double value = tangent * A;
	//	double vas1 = -y * value;
	//	double vas2 = z * value;
	//	double vas1as2 = vas1 * z;

	//	kInitialData[0] += value;
	//	kInitialData[1] += vas1;
	//	kInitialData[2] += vas2;

	//	kInitialData[5] += vas1 * -y;
	//	kInitialData[6] += vas1as2;

	//	kInitialData[10] += vas2 * z;
	//}

	//kInitialData[4] = kInitialData[1];
	//kInitialData[8] = kInitialData[2];
	//kInitialData[9] = kInitialData[6];


	return kInitial;
}

const Vector&
FiberSectionShearWall::getSectionDeformation(void)
{
	return e;
}

const Matrix&
FiberSectionShearWall::getSectionTangent(void)
{
	return *ks;
}

const Vector&
FiberSectionShearWall::getStressResultant(void)
{
	return *s;
}

SectionForceDeformation*
FiberSectionShearWall::getCopy(void)
{
	FiberSectionShearWall *theCopy = new FiberSectionShearWall();
	theCopy->setTag(this->getTag());

	theCopy->numFibers = numFibers;
	theCopy->sizeFibers = numFibers;

	if (numFibers != 0) {
		theCopy->theMaterials = new NDMaterial *[numFibers];

		if (theCopy->theMaterials == 0) {
			opserr << "FiberSectionShearWall::FiberSectionShearWall -- failed to allocate Material pointers\n";
			exit(-1);
		}

		theCopy->matData = new double[numFibers * 3];

		if (theCopy->matData == 0) {
			opserr << "FiberSectionShearWall::FiberSectionShearWall -- failed to allocate double array for material data\n";
			exit(-1);
		}


		for (int i = 0; i < numFibers; i++) {
			theCopy->matData[i * 3] = matData[i * 3];
			theCopy->matData[i * 3 + 1] = matData[i * 3 + 1];
			theCopy->matData[i * 3 + 2] = matData[i * 3 + 2];
			theCopy->theMaterials[i] = theMaterials[i]->getCopy();

			if (theCopy->theMaterials[i] == 0) {
				opserr << "FiberSectionShearWall::getCopy -- failed to get copy of a Material\n";
				exit(-1);
			}
		}
	}

	theCopy->e = e;
	theCopy->QzBar = QzBar;
	theCopy->QyBar = QyBar;
	theCopy->Abar = Abar;
	theCopy->yBar = yBar;
	theCopy->zBar = zBar;

	for (int i = 0; i < 64; i++)
		theCopy->kData[i] = kData[i];

	for (int i = 0; i < 8; i++)
		theCopy->sData[i] = kData[i];

	if (sectionIntegr != 0)
		theCopy->sectionIntegr = sectionIntegr->getCopy();
	else
		theCopy->sectionIntegr = 0;

	return theCopy;
}

const ID&
FiberSectionShearWall::getType()
{
	return code;
}

int
FiberSectionShearWall::getOrder() const
{
	return 4;
}

int
FiberSectionShearWall::commitState(void)
{
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->commitState();

	return err;
}

int
FiberSectionShearWall::revertToLastCommit(void)
{
	int err = 0;

	s->Zero();
	ks->Zero();

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
			zLocs[i] = matData[3 * i + 1];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double z = zLocs[i] - zBar;
		double A = fiberArea[i];

		// invoke revertToLast on the material
		err += theMat->revertToLastCommit();

		const Vector &stress = theMat->getStress();
		const Matrix &tangent = theMat->getTangent();

		Matrix crdCoeff(3, 8);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(0, 2) = -z;
		crdCoeff(0, 3) = -y * z;
		crdCoeff(1, 4) = 1;
		crdCoeff(1, 5) = -z;
		crdCoeff(2, 6) = 1;
		crdCoeff(2, 7) = -y;

		Matrix secStiff(8, 8), CoeffD(3, 8);
		Vector secForce(8);

		secForce.addMatrixTransposeVector(1.0, crdCoeff, stress, 1.0);
		CoeffD.addMatrixProduct(0.0, tangent, crdCoeff, 1.0);
		secStiff.addMatrixTransposeProduct(0.0, crdCoeff, CoeffD, 1.0);
		(*s) += secForce * A;
		(*ks) += secStiff * A;
	}
	return err;
}

int
FiberSectionShearWall::revertToStart(void)
{
	// revert the fibers to start    
	int err = 0;

	s->Zero();
	ks->Zero();

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
			zLocs[i] = matData[3 * i + 1];
			fiberArea[i] = matData[3 * i + 2];
		}
	}

	for (int i = 0; i < numFibers; i++) {
		NDMaterial *theMat = theMaterials[i];
		double y = yLocs[i] - yBar;
		double z = zLocs[i] - zBar;
		double A = fiberArea[i];

		// invoke revertToStart on the material
		// invoke revertToLast on the material
		err += theMat->revertToStart();

		const Vector &stress = theMat->getStress();
		const Matrix &tangent = theMat->getTangent();

		Matrix crdCoeff(3, 8);
		crdCoeff(0, 0) = 1;
		crdCoeff(0, 1) = y;
		crdCoeff(0, 2) = -z;
		crdCoeff(0, 3) = -y * z;
		crdCoeff(1, 4) = 1;
		crdCoeff(1, 5) = -z;
		crdCoeff(2, 6) = 1;
		crdCoeff(2, 7) = -y;

		Matrix secStiff(8, 8), CoeffD(3, 8);
		Vector secForce(8);

		secForce.addMatrixTransposeVector(1.0, crdCoeff, stress, 1.0);
		CoeffD.addMatrixProduct(0.0, tangent, crdCoeff, 1.0);
		secStiff.addMatrixTransposeProduct(0.0, crdCoeff, CoeffD, 1.0);
		(*s) += secForce * A;
		(*ks) += secStiff * A;
	}

	return err;
}

int
FiberSectionShearWall::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	//// create an id to send objects tag and numFibers, 
	////     size 3 so no conflict with matData below if just 1 fiber
	//static ID data(3);
	//data(0) = this->getTag();
	//data(1) = numFibers;
	//int dbTag = this->getDbTag();
	//theTorsion->setDbTag(dbTag);
	//data(2) = theTorsion->getClassTag();

	//res += theChannel.sendID(dbTag, commitTag, data);
	//if (res < 0) {
	//	opserr << "FiberSectionShearWall::sendSelf - failed to send ID data\n";
	//	return res;
	//}

	//theTorsion->sendSelf(commitTag, theChannel);

	//if (numFibers != 0) {

	//	// create an id containingg classTag and dbTag for each material & send it
	//	ID materialData(2 * numFibers);
	//	for (int i = 0; i < numFibers; i++) {
	//		UniaxialMaterial *theMat = theMaterials[i];
	//		materialData(2 * i) = theMat->getClassTag();
	//		int matDbTag = theMat->getDbTag();
	//		if (matDbTag == 0) {
	//			matDbTag = theChannel.getDbTag();
	//			if (matDbTag != 0)
	//				theMat->setDbTag(matDbTag);
	//		}
	//		materialData(2 * i + 1) = matDbTag;
	//	}

	//	res += theChannel.sendID(dbTag, commitTag, materialData);
	//	if (res < 0) {
	//		opserr << "FiberSectionShearWall::sendSelf - failed to send material data\n";
	//		return res;
	//	}

	//	// send the fiber data, i.e. area and loc
	//	Vector fiberData(matData, 3 * numFibers);
	//	res += theChannel.sendVector(dbTag, commitTag, fiberData);
	//	if (res < 0) {
	//		opserr << "FiberSectionShearWall::sendSelf - failed to send fiber data\n";
	//		return res;
	//	}

	//	// now invoke send(0 on all the materials
	//	for (int j = 0; j < numFibers; j++)
	//		theMaterials[j]->sendSelf(commitTag, theChannel);
	//}

	return res;
}

int
FiberSectionShearWall::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	//static ID data(3);

	//int dbTag = this->getDbTag();
	//res += theChannel.recvID(dbTag, commitTag, data);

	//if (res < 0) {
	//	opserr << "FiberSectionShearWall::recvSelf - failed to recv ID data\n";
	//	return res;
	//}

	//this->setTag(data(0));

	//if (theTorsion == 0) {
	//	int cTag = data(2);
	//	theTorsion = theBroker.getNewUniaxialMaterial(cTag);
	//	theTorsion->setDbTag(dbTag);
	//}
	//if (theTorsion == 0) {
	//	opserr << "FiberSectionShearWall::recvSelf - failed to get torsion material \n";
	//	return -1;
	//}
	//if (theTorsion->recvSelf(commitTag, theChannel, theBroker) < 0) {
	//	opserr << "FiberSectionShearWall::recvSelf - torsion failed to recvSelf \n";
	//	return -2;
	//}

	//// recv data about materials objects, classTag and dbTag
	//if (data(1) != 0) {
	//	ID materialData(2 * data(1));
	//	res += theChannel.recvID(dbTag, commitTag, materialData);
	//	if (res < 0) {
	//		opserr << "FiberSectionShearWall::recvSelf - failed to recv material data\n";
	//		return res;
	//	}

	//	// if current arrays not of correct size, release old and resize
	//	if (theMaterials == 0 || numFibers != data(1)) {
	//		// delete old stuff if outa date
	//		if (theMaterials != 0) {
	//			for (int i = 0; i < numFibers; i++)
	//				delete theMaterials[i];
	//			delete[] theMaterials;
	//			if (matData != 0)
	//				delete[] matData;
	//			matData = 0;
	//			theMaterials = 0;
	//		}

	//		// create memory to hold material pointers and fiber data
	//		numFibers = data(1);
	//		sizeFibers = data(1);
	//		if (numFibers != 0) {

	//			theMaterials = new UniaxialMaterial *[numFibers];

	//			if (theMaterials == 0) {
	//				opserr << "FiberSectionShearWall::recvSelf -- failed to allocate Material pointers\n";
	//				exit(-1);
	//			}

	//			for (int j = 0; j < numFibers; j++)
	//				theMaterials[j] = 0;

	//			matData = new double[numFibers * 3];

	//			if (matData == 0) {
	//				opserr << "FiberSectionShearWall::recvSelf  -- failed to allocate double array for material data\n";
	//				exit(-1);
	//			}
	//		}
	//	}

	//	Vector fiberData(matData, 3 * numFibers);
	//	res += theChannel.recvVector(dbTag, commitTag, fiberData);
	//	if (res < 0) {
	//		opserr << "FiberSectionShearWall::recvSelf - failed to recv fiber data\n";
	//		return res;
	//	}

	//	int i;
	//	for (i = 0; i < numFibers; i++) {
	//		int classTag = materialData(2 * i);
	//		int dbTag = materialData(2 * i + 1);

	//		// if material pointed to is blank or not of corrcet type, 
	//		// release old and create a new one
	//		if (theMaterials[i] == 0)
	//			theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);
	//		else if (theMaterials[i]->getClassTag() != classTag) {
	//			delete theMaterials[i];
	//			theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);
	//		}

	//		if (theMaterials[i] == 0) {
	//			opserr << "FiberSectionShearWall::recvSelf -- failed to allocate double array for material data\n";
	//			exit(-1);
	//		}

	//		theMaterials[i]->setDbTag(dbTag);
	//		res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
	//	}

	//	QzBar = 0.0;
	//	QyBar = 0.0;
	//	Abar = 0.0;
	//	double yLoc, zLoc, Area;

	//	// Recompute centroid
	//	for (i = 0; i < numFibers; i++) {
	//		yLoc = matData[3 * i];
	//		zLoc = matData[3 * i + 1];
	//		Area = matData[3 * i + 2];
	//		Abar += Area;
	//		QzBar += yLoc * Area;
	//		QyBar += zLoc * Area;
	//	}

	//	yBar = QzBar / Abar;
	//	zBar = QyBar / Abar;
	//}

	return res;
}

void
FiberSectionShearWall::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "\tSection type: FiberSectionShearWall, tag: " << this->getTag() << endln;
		s << "\tNumber of Fibers: " << numFibers << endln;
		s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;

		s << "\nMaterial Information :";
		if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
			for (int i = 0; i < numFibers; i++) {
				s << "\n\tFiber  Tag: " << i + 1;
				s << "\n\tLocation (y, z) = (" << matData[3 * i] << ", " << matData[3 * i + 1] << ")";
				s << "\n\tArea = " << matData[3 * i + 2] << endln;
				theMaterials[i]->Print(s, flag);

			}
		}
	}
}

Response*
FiberSectionShearWall::setResponse(const char **argv, int argc, OPS_Stream &output)
{

	//const ID &type = this->getType();
	//int typeSize = this->getOrder();

	Response *theResponse = 0;


	if (argc > 2 || strcmp(argv[0], "fiber") == 0) {

		int key = numFibers;
		int passarg = 2;

		if (argc == 4) {  // fiber number was input directly
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
		}

		if (key < numFibers && key >= 0) {

			theResponse = theMaterials[key]->setResponse(&argv[passarg], argc - passarg, output);
	
		}
	}
	else
	{
		opserr << "this recorder type is not supported in SW-fiber section " << endln;
	}

	return theResponse;
}


int
FiberSectionShearWall::getResponse(int responseID, Information &sectInfo)
{
	//// Just call the base class method ... don't need to define
	//// this function, but keeping it here just for clarity
	//if (responseID == 5) {
	//	int numData = 5 * numFibers;
	//	Vector data(numData);
	//	int count = 0;
	//	for (int j = 0; j < numFibers; j++) {
	//		double yLoc, zLoc, A, stress, strain;
	//		yLoc = matData[3 * j];
	//		zLoc = matData[3 * j + 1];
	//		A = matData[3 * j + 2];
	//		stress = theMaterials[j]->getStress();
	//		strain = theMaterials[j]->getStrain();
	//		data(count) = yLoc; data(count + 1) = zLoc; data(count + 2) = A;
	//		data(count + 3) = stress; data(count + 4) = strain;
	//		count += 5;
	//	}
	//	return sectInfo.setVector(data);
	//}
	//else  if (responseID == 6) {
	//	int count = 0;
	//	for (int j = 0; j < numFibers; j++) {
	//		if (theMaterials[j]->hasFailed() == true)
	//			count++;
	//	}
	//	return sectInfo.setInt(count);
	//}
	//else  if (responseID == 7) {
	//	int count = 0;
	//	for (int j = 0; j < numFibers; j++) {
	//		if (theMaterials[j]->hasFailed() == true) {
	//			count += 1;
	//		}
	//	}
	//	if (count == numFibers)
	//		count = 1;
	//	else
	//		count = 0;

	//	return sectInfo.setInt(count);
	//}
	//else  if (responseID == 10) {

	//	return sectInfo.setDouble(getEnergy());
	//}

	//return SectionForceDeformation::getResponse(responseID, sectInfo);
	return 0;
}

int
FiberSectionShearWall::setParameter(const char **argv, int argc, Parameter &param)
{
	//if (argc < 1)
	//	return -1;

	int result = 0;

	//// A material parameter
	//if (strstr(argv[0], "material") != 0) {

	//	// Get the tag of the material
	//	int paramMatTag = atoi(argv[1]);

	//	// Loop over fibers to find the right material(s)
	//	int ok = 0;
	//	for (int i = 0; i < numFibers; i++)
	//		if (paramMatTag == theMaterials[i]->getTag()) {
	//			ok = theMaterials[i]->setParameter(&argv[2], argc - 2, param);
	//			if (ok != -1)
	//				result = ok;
	//		}

	//	if (paramMatTag == theTorsion->getTag()) {
	//		ok = theTorsion->setParameter(&argv[2], argc - 2, param);
	//		if (ok != -1)
	//			result = ok;
	//	}
	//	return result;
	//}

	//// Check if it belongs to the section integration
	//else if (strstr(argv[0], "integration") != 0) {
	//	if (sectionIntegr != 0)
	//		return sectionIntegr->setParameter(&argv[1], argc - 1, param);
	//	else
	//		return -1;
	//}

	//int ok = 0;

	//// loop over every material
	//for (int i = 0; i < numFibers; i++) {
	//	ok = theMaterials[i]->setParameter(argv, argc, param);
	//	if (ok != -1)
	//		result = ok;
	//}

	//// Don't really need to do this in "default" mode
	////ok = theTorsion->setParameter(argv, argc, param);
	////if (ok != -1)
	////  result = ok;

	//if (sectionIntegr != 0) {
	//	ok = sectionIntegr->setParameter(argv, argc, param);
	//	if (ok != -1)
	//		result = ok;
	//}

	return result;
}
