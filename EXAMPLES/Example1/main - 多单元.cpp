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

// $Revision: 1.6 $
// $Date: 2009-11-02 22:23:58 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/Example1/main.cpp,v $


// File: ~/model/main.C
//
// Written: fmk 08/99
//
// Purpose: this file contains a C++ main procedure to perform the analysis
// of example1 (found in most documents). In the main() procedure:
// 	1) each object of the domain, i.e. Nodes, Elements, Constraints,
//	   and LoadPattern objects are created and then added to the Domain.
//	2) the components of the analysis object are constructed and then
//	   the Analysis object is created.
//	3) the analysis is performed.
//	4) the results are printed - here the contents of Domain and end of
//	   the analysis operation.

// standard C++ includes

#include <stdlib.h>

#include <OPS_Globals.h>
#include <StandardStream.h>

#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <Node.h>
#include <Truss.h>
#include <ElasticMaterial.h>
#include <SP_Constraint.h>
#include <LoadPattern.h>
#include <LinearSeries.h>
#include <NodalLoad.h>

// includes for the analysis classes
#include <StaticAnalysis.h>
#include <AnalysisModel.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <DOF_Numberer.h>
#include <RCM.h>
#include <LoadControl.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <CTestNormDispIncr.h>

#include <ovalConcrete.h>
#include <steelTimo.h>
//#include <NDFiber.h>
#include <NDFiber2D.h>
#include <FiberSectionTimo.h>
#include <TimoshenkoModified.h>

#include <math.h>

#include <string>

// init the global variabled defined in OPS_Globals.h
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;


int main(int argc, char **argv)
{

	int analysisStep = 0;

	Domain *theDomain = new Domain();

	//////////// build model

	/////////////////  beam
	Node *node1 = new Node(1, 3, 0.0, 0.0);
	Node *node2 = new Node(2, 3, 0.0, 500.0);
	Node *node3 = new Node(3, 3, 0.0, 1000.0);
	Node *node4 = new Node(4, 3, 0.0, 1500.0);
	Node *node8 = new Node(8, 3, 0.0, 2000.0);

	theDomain->addNode(node1);
	theDomain->addNode(node2);
	theDomain->addNode(node3);
	theDomain->addNode(node4);
	theDomain->addNode(node8);

	SP_Constraint *sp1 = new SP_Constraint(1, 0, 0.0, true);
	SP_Constraint *sp2 = new SP_Constraint(1, 1, 0.0, true);
	SP_Constraint *sp3 = new SP_Constraint(1, 2, 0.0, true);
	theDomain->addSP_Constraint(sp1);
	theDomain->addSP_Constraint(sp2);
	theDomain->addSP_Constraint(sp3);

	NDMaterial *concrete = new ovalConcrete(1, -40.0, -2.0e-3, -2.0e-2, 1.0e-4, 1.0e-3, 0.9, 0.7);
	NDMaterial *steel = new steelTimo(2, 336.0, 2.0e5, 0.035, 18, 0.925, 0.15);
	UniaxialMaterial *elastic = new ElasticMaterial(3, 1.0e12);

	Fiber **fibers = new Fiber*[14];
	fibers[0]  = new NDFiber2d(1,  *concrete, 12622.5, -453.75);
	fibers[1]  = new NDFiber2d(2,  *concrete, 12622.5, -371.25);
	fibers[2]  = new NDFiber2d(3,  *concrete, 12622.5, -288.75);
	fibers[3]  = new NDFiber2d(4,  *concrete, 12622.5, -206.25);
	fibers[4]  = new NDFiber2d(5,  *concrete, 12622.5, -123.75);
	fibers[5]  = new NDFiber2d(6,  *concrete, 12622.5, -41.250);
	fibers[6]  = new NDFiber2d(7,  *concrete, 12622.5, 41.250);
	fibers[7]  = new NDFiber2d(8,  *concrete, 12622.5, 123.75);
	fibers[8]  = new NDFiber2d(9,  *concrete, 12622.5, 206.25);
	fibers[9]  = new NDFiber2d(10, *concrete, 12622.5, 288.75);
	fibers[10] = new NDFiber2d(11, *concrete, 12622.5, 371.25);
	fibers[11] = new NDFiber2d(12, *concrete, 12622.5, 453.75);
	fibers[12] = new NDFiber2d(13, *steel, 1300.0, -440.0);
	fibers[13] = new NDFiber2d(14, *steel, 1300.0, 440.0);

	SectionForceDeformation* theSection = new FiberSectionTimo(1, 14, fibers);

	TimoshenkoModified* beam1 = new TimoshenkoModified(1, 1, 2, theSection, 0.0);
	TimoshenkoModified* beam2 = new TimoshenkoModified(2, 2, 3, theSection, 0.0);
	TimoshenkoModified* beam3 = new TimoshenkoModified(3, 3, 4, theSection, 0.0);
	TimoshenkoModified* beam4 = new TimoshenkoModified(4, 4, 8, theSection, 0.0);

	theDomain->addElement(beam1);
	theDomain->addElement(beam2);
	theDomain->addElement(beam3);
	theDomain->addElement(beam4);

	/////////////////  truss
//	Node *node12 = new Node(12, 3, 0.0, 2001.0);
	Node *node22 = new Node(22, 3, 1.0, 2000.0);

//	theDomain->addNode(node12);
	theDomain->addNode(node22);

	//SP_Constraint *sp4 = new SP_Constraint(12, 0, 0.0, true);
	//SP_Constraint *sp5 = new SP_Constraint(12, 1, 0.0, true);
	//SP_Constraint *sp6 = new SP_Constraint(12, 2, 0.0, true);
	SP_Constraint *sp7 = new SP_Constraint(22, 0, 0.0, true);
	SP_Constraint *sp8 = new SP_Constraint(22, 1, 0.0, true);
	SP_Constraint *sp9 = new SP_Constraint(22, 2, 0.0, true);

	//theDomain->addSP_Constraint(sp4);
	//theDomain->addSP_Constraint(sp5);
	//theDomain->addSP_Constraint(sp6);
	theDomain->addSP_Constraint(sp7);
	theDomain->addSP_Constraint(sp8);
	theDomain->addSP_Constraint(sp9);

//	Truss *truss1 = new Truss(12, 2, 8, 12, *elastic, 1.0);
	Truss *truss2 = new Truss(22, 2, 8, 22, *elastic, 1.0);

//	theDomain->addElement(truss1);
	theDomain->addElement(truss2);

	//////////// add vertical loads

	TimeSeries *theSeries1 = new LinearSeries();

	LoadPattern *theLoadPattern1 = new LoadPattern(1);

	theLoadPattern1->setTimeSeries(theSeries1);
	theDomain->addLoadPattern(theLoadPattern1);

	Vector theLoadValues1(3);
	theLoadValues1(1) = 1.0;
	NodalLoad *theLoad1 = new NodalLoad(1, 8, theLoadValues1);
	theDomain->addNodalLoad(theLoad1, 1);

	//////////////  gravity analysis

	double incrFacotr1 = -7.0e3;
	AnalysisModel     *theModel1 = new AnalysisModel();
	EquiSolnAlgo      *theSolnAlgo1 = new NewtonRaphson();
	StaticIntegrator  *theIntegrator1 = new LoadControl(incrFacotr1, 1, incrFacotr1, incrFacotr1);
	ConstraintHandler *theHandler1 = new TransformationConstraintHandler();
	RCM               *theRCM1 = new RCM();
	DOF_Numberer      *theNumberer1 = new DOF_Numberer(*theRCM1);
	BandGenLinSolver  *theSolver1 = new BandGenLinLapackSolver();
	LinearSOE         *theSOE1 = new BandGenLinSOE(*theSolver1);
	ConvergenceTest   *theTest1 = new CTestNormDispIncr(1.0e-8, 30, 0);

	StaticAnalysis    theAnalysis1(*theDomain,
		*theHandler1,
		*theNumberer1,
		*theModel1,
		*theSolnAlgo1,
		*theSOE1,
		*theIntegrator1,
		theTest1);

	theAnalysis1.analyze(100);
	analysisStep += 100;

	opserr << node8->getDisp()[0] << "\t" << node8->getDisp()[1] << "\t" << node8->getDisp()[2] << endln;
	///// const load

	theDomain->setLoadConstant();
	theDomain->setCurrentTime(0.0);
	theDomain->setCommittedTime(0.0);

	////////////// add push over loads

	TimeSeries *theSeries2 = new LinearSeries();

	LoadPattern *theLoadPattern2 = new LoadPattern(2);

	theLoadPattern2->setTimeSeries(theSeries2);
	theDomain->addLoadPattern(theLoadPattern2);

	Vector theLoadValues2(3);
	theLoadValues2(0) = 1.0e12;
	NodalLoad *theLoad2 = new NodalLoad(2, 8, theLoadValues2);

	theDomain->addNodalLoad(theLoad2, 2);

	//////////////  push over analysis

	double incrFacotr2 = 1.0e-1;

	int totalIterNum = 1;
	Vector maxIterDisp(totalIterNum);
	maxIterDisp(0) = 71.28;

	ofstream file;
	//file.open("trailStrainStress.txt");
	file.open("info.txt");
	for (int i = 0; i < totalIterNum; i++)
	{
		double theMaxIterDisp = maxIterDisp[i];
		int maxStep = 4 * floor(abs(theMaxIterDisp / incrFacotr2));
		int count = 0;
		while (count < maxStep)
		{
			int multiFactor;
			if ((count < maxStep / 4) || (count < maxStep && count >= maxStep * 3 / 4))
				multiFactor = 1.0;
			else
				multiFactor = -1.0;

			incrFacotr2 = abs(incrFacotr2) * multiFactor;

			AnalysisModel     *theModel2 = new AnalysisModel();
			EquiSolnAlgo      *theSolnAlgo2 = new NewtonRaphson();
			StaticIntegrator  *theIntegrator2 = new LoadControl(incrFacotr2, 1, incrFacotr2, incrFacotr2);
			ConstraintHandler *theHandler2 = new TransformationConstraintHandler();
			RCM               *theRCM2 = new RCM();
			DOF_Numberer      *theNumberer2 = new DOF_Numberer(*theRCM2);
			BandGenLinSolver  *theSolver2 = new BandGenLinLapackSolver();
			LinearSOE         *theSOE2 = new BandGenLinSOE(*theSolver2);
			ConvergenceTest   *theTest2 = new CTestNormDispIncr(1.0e-8, 30, 0);

			StaticAnalysis    theAnalysis2(*theDomain,
				*theHandler2,
				*theNumberer2,
				*theModel2,
				*theSolnAlgo2,
				*theSOE2,
				*theIntegrator2,
				theTest2);

			int condition = theAnalysis2.analyze(1);
			if (condition != 0)
				break;
			count++;
			analysisStep++;

			opserr << analysisStep << endln;
		}
	}
	file.close();


	///// perturb

	int halfPertStep = 100;


	/////////////////////////////////

	int numEle = 1;
	int numFibers = 5;

	FiberSectionTimo **theSec = new FiberSectionTimo *[numEle];
	theSec[0] = (FiberSectionTimo *)beam1->theSection;
	theSec[1] = (FiberSectionTimo *)beam2->theSection;
	theSec[2] = (FiberSectionTimo *)beam3->theSection;
	theSec[3] = (FiberSectionTimo *)beam4->theSection;

	ofstream file1, file2, file3;

	file1.precision(12);
	file2.precision(12);
	file3.precision(12);

	//for (int count = 0; count < 4; count++)
	//{

	//	std::string s1;

	//	s1 = "D:\\pertDisp" + std::to_string(count) + ".in";

	//	file1.open(s1);

	//	double pertVertDispValue = 1.0e-10;
	//	double pertHorzDispValue = 1.0e-10;
	//	double pertRotateValue = 1.0e-12;

	//	Vector secStrain = theSec[count]->getSectionDeformation();

	//	Vector secStrainTr(3), secForce;
	//	for (int i = 0; i < 2 * halfPertStep + 1; i++)
	//	{
	//		secStrainTr(0) = secStrain(0) - halfPertStep * pertVertDispValue + i * pertVertDispValue;
	//		secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + i * pertRotateValue;
	//		secStrainTr(2) = secStrain(2) - halfPertStep * pertHorzDispValue + i * pertHorzDispValue;

	//		file1 << secStrainTr(0) << "\t" << secStrainTr(1) << "\t" << secStrainTr(2) << "\t" << endln;
	//	}
	//	file1.close();
	//}


	//for (int count = 0; count < 4; count++)
	//{
	//	std::string s1, s2, s3;

	//	s1 = "D:\\normFPertDvDh" + std::to_string(count) + ".out";
	//	s2 = "D:\\momentPertDvDh" + std::to_string(count) + ".out";
	//	s3 = "D:\\shearFPertDvDh" + std::to_string(count) + ".out";

	//	file1.open(s1);
	//	file2.open(s2);
	//	file3.open(s3);

	//	double pertVertDispValue = 1.0e-10;
	//	double pertHorzDispValue = 1.0e-10;

	//	Vector secStrain = theSec[count]->getSectionDeformation();

	//	Vector secStrainTr(3), secForce;

	//	for (int i = 0; i < 2 * halfPertStep + 1; i++)
	//	{
	//		for (int j = 0; j < 2 * halfPertStep + 1; j++)
	//		{
	//				secStrainTr(0) = secStrain(0) - halfPertStep * pertVertDispValue + i * pertVertDispValue;
	//				secStrainTr(2) = secStrain(2) - halfPertStep * pertHorzDispValue + j * pertHorzDispValue;

	//				theSec[count]->setTrialSectionDeformation(secStrainTr);

	//				secForce = theSec[count]->getStressResultant();

	//				file1 << secForce(0) << "\t";
	//				file2 << secForce(1) << "\t";
	//				file3 << secForce(2) << "\t";
	//		}
	//		file1 << endln;
	//		file2 << endln;
	//		file3 << endln;

	//	}

	//	file1.close();
	//	file2.close();
	//	file3.close();
	//}


	//for (int count = 0; count < 4; count++)
	//{
	//	std::string s1, s2, s3;

	//	s1 = "D:\\normFPertDvDt" + std::to_string(count) + ".out";
	//	s2 = "D:\\momentPertDvDt" + std::to_string(count) + ".out";
	//	s3 = "D:\\shearFPertDvDt" + std::to_string(count) + ".out";

	//	file1.open(s1);
	//	file2.open(s2);
	//	file3.open(s3);

	//	double pertVertDispValue = 1.0e-10;
	//	double pertRotateValue = 1.0e-12;

	//	Vector secStrain = theSec[count]->getSectionDeformation();

	//	Vector secStrainTr(3), secForce;

	//	for (int i = 0; i < 2 * halfPertStep + 1; i++)
	//	{
	//		for (int j = 0; j < 2 * halfPertStep + 1; j++)
	//		{
	//			secStrainTr(0) = secStrain(0) - halfPertStep * pertVertDispValue + i * pertVertDispValue;
	//			secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + j * pertRotateValue;

	//			theSec[count]->setTrialSectionDeformation(secStrainTr);

	//			secForce = theSec[count]->getStressResultant();

	//			file1 << secForce(0) << "\t";
	//			file2 << secForce(1) << "\t";
	//			file3 << secForce(2) << "\t";
	//		}
	//		file1 << endln;
	//		file2 << endln;
	//		file3 << endln;

	//	}
	//	file1.close();
	//	file2.close();
	//	file3.close();
	//}


	//for (int count = 0; count < 4; count++)
	//{
	//	std::string s1, s2, s3;

	//	s1 = "D:\\normFPertDhDt" + std::to_string(count) + ".out";
	//	s2 = "D:\\momentPertDhDt" + std::to_string(count) + ".out";
	//	s3 = "D:\\shearFPertDhDt" + std::to_string(count) + ".out";

	//	file1.open(s1);
	//	file2.open(s2);
	//	file3.open(s3);

	//	double pertHorzispValue = 1.0e-10;
	//	double pertRotateValue = 1.0e-12;

	//	Vector secStrain = theSec[count]->getSectionDeformation();

	//	Vector secStrainTr(3), secForce;

	//	for (int i = 0; i < 2 * halfPertStep + 1; i++)
	//	{
	//		for (int j = 0; j < 2 * halfPertStep + 1; j++)
	//		{
	//			secStrainTr(2) = secStrain(2) - halfPertStep * pertHorzispValue + i * pertHorzispValue;
	//			secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + j * pertRotateValue;

	//			theSec[count]->setTrialSectionDeformation(secStrainTr);

	//			secForce = theSec[count]->getStressResultant();

	//			file1 << secForce(0) << "\t";
	//			file2 << secForce(1) << "\t";
	//			file3 << secForce(2) << "\t";
	//		}
	//		file1 << endln;
	//		file2 << endln;
	//		file3 << endln;

	//	}

	//	file1.close();
	//	file2.close();
	//	file3.close();
	//}

	///////////////////////////////////



	//NDMaterial **matPtr = new NDMaterial *[48];
	//for (int i = 0; i < numEle; i++)
	//	for (int j = 0; j < numFibers; j++)
	//		matPtr[i * numFibers + j] = theSec[i]->theMaterials[j];

	//for (int count = 0; count < numFibers * numEle; count++)
	//{
	//	std::string s1, s2, s3, s4;

	//	s1 = "D:\\pertNormStress"  + std::to_string(count) + ".out";
	//	s2 = "D:\\pertShearStress" + std::to_string(count) + ".out";
	//	s3 = "D:\\pertNormStrain"  + std::to_string(count) + ".in";
	//	s4 = "D:\\pertShearStrain" + std::to_string(count) + ".in";

	//	file1.open(s1);
	//	file2.open(s2);
	//	file3.open(s3);
	//	file4.open(s4);

	//	double pertNormStrain = 1.0e-9;
	//	double pertShearStrain = 1.0e-8;

	//	Vector strain = matPtr[count]->getStrain();
	//	for (int i = 0; i < 2 * halfPertStep + 1; i++)
	//	{
	//		Vector strainTr(2), stress;

	//		for (int j = 0; j < 2 * halfPertStep + 1; j++)
	//		{
	//			strainTr(0) = strain(0) - halfPertStep * pertNormStrain + i * pertNormStrain;
	//			strainTr(1) = strain(1) - halfPertStep * pertShearStrain + j * pertShearStrain;

	//			matPtr[count]->setTrialStrain(strainTr);

	//			stress = matPtr[count]->getStress();

	//			file1 << stress(0) << "\t";
	//			file2 << stress(1) << "\t";
	//			if (i == 0)
	//				file4 << strainTr(1) << endln;
	//		}
	//		file1 << endln;
	//		file2 << endln;
	//		file3 << strainTr(0) << endln;
	//	}
	//	file1.close();
	//	file2.close();
	//	file3.close();
	//	file4.close();
	//}


	NDMaterial **matPtr = new NDMaterial *[5];
	for (int j = 0; j < numFibers; j++)
		matPtr[j] = theSec[0]->theMaterials[j];

	for (int count = 0; count < 12; count++)
	{
		std::string s1;

		s1 = "D:\\pertStrain" + std::to_string(count) + ".in";

		file1.open(s1);

		Vector strain = matPtr[count]->getStrain();

		Vector strainTr(2);

		double pertNormStrain = 1.0e-9;
		double pertShearStrain = 1.0e-8;
		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			strainTr(0) = strain(0) - halfPertStep * pertNormStrain + i * pertNormStrain;
			strainTr(1) = strain(1) - halfPertStep * pertShearStrain + i * pertShearStrain;
			file1 << strainTr(0) << "\t" << strainTr(1) << endln;
		}
		file1.close();
	}

	for (int count = 0; count < 12; count++)
	{
		std::string s1, s2, s3, s4;

		s1 = "D:\\normStressPert" + std::to_string(count) + ".out";
		s2 = "D:\\shearStressPert" + std::to_string(count) + ".out";

		file1.open(s1);
		file2.open(s2);

		double pertNormStrain = 1.0e-8;
		double pertShearStrain = 1.0e-7;

		Vector strain = matPtr[count]->getStrain();
		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			Vector strainTr(2), stress;

			for (int j = 0; j < 2 * halfPertStep + 1; j++)
			{
				strainTr(0) = strain(0) - halfPertStep * pertNormStrain + i * pertNormStrain;
				strainTr(1) = strain(1) - halfPertStep * pertShearStrain + j * pertShearStrain;

				matPtr[count]->setTrialStrain(strainTr);

				stress = matPtr[count]->getStress();

				file1 << stress(0) << "\t";
				file2 << stress(1) << "\t";
			}
			file1 << endln;
			file2 << endln;
		}
		file1.close();
		file2.close();
	}



	opserr << "finished." << endln;
	getchar();

	return 0;
}

