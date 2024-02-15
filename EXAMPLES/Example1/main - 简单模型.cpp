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
	Node *node2 = new Node(2, 3, 0.0, 100.0);
	//Node *node3 = new Node(3, 3, 0.0, 1000.0);
	//Node *node4 = new Node(4, 3, 0.0, 1500.0);
	Node *node8 = new Node(8, 3, 0.0, 2000.0);

	theDomain->addNode(node1);
	theDomain->addNode(node2);
	//theDomain->addNode(node3);
	//theDomain->addNode(node4);
	theDomain->addNode(node8);

	SP_Constraint *sp1 = new SP_Constraint(1, 0, 0.0, true);
	SP_Constraint *sp2 = new SP_Constraint(1, 1, 0.0, true);
	SP_Constraint *sp3 = new SP_Constraint(1, 2, 0.0, true);
	theDomain->addSP_Constraint(sp1);
	theDomain->addSP_Constraint(sp2);
	theDomain->addSP_Constraint(sp3);

	NDMaterial *concrete = new ovalConcrete(1, -40.0, -2.0e-3, -2.0e-2, 1.0e-4, 1.0e-3, 0.9, 0.2);
	NDMaterial *steel = new steelTimo(2, 336.0, 2.0e5, 0.035, 18, 0.925, 0.15);
	UniaxialMaterial *elastic = new ElasticMaterial(3, 1.0e12);

	Fiber **fibers = new Fiber*[6];
	fibers[0] = new NDFiber2d(1, *concrete, 37867.5, -371.25);
	fibers[1] = new NDFiber2d(2, *concrete, 37867.5, -123.75);
	fibers[2] = new NDFiber2d(3, *concrete, 37867.5, 123.75);
	fibers[3] = new NDFiber2d(4, *concrete, 37867.5, 371.25);
	fibers[4] = new NDFiber2d(5, *steel, 1300.0, -440.0);
	fibers[5] = new NDFiber2d(6, *steel, 1300.0, 440.0);

	SectionForceDeformation* theSection = new FiberSectionTimo(1, 6, fibers);

	TimoshenkoModified* beam1 = new TimoshenkoModified(1, 1, 2, theSection, 0.0);
	//TimoshenkoModified* beam2 = new TimoshenkoModified(2, 2, 3, theSection, 0.0);
	//TimoshenkoModified* beam3 = new TimoshenkoModified(3, 3, 4, theSection, 0.0);
	TimoshenkoModified* beam4 = new TimoshenkoModified(4, 2, 8, theSection, 0.0);

	theDomain->addElement(beam1);
	//theDomain->addElement(beam2);
	//theDomain->addElement(beam3);
	theDomain->addElement(beam4);

	/////////////////  truss
	Node *node12 = new Node(12, 3, 0.0, 2001.0);
	Node *node22 = new Node(22, 3, 1.0, 2000.0);

	theDomain->addNode(node12);
	theDomain->addNode(node22);

	SP_Constraint *sp4 = new SP_Constraint(12, 0, 0.0, true);
	SP_Constraint *sp5 = new SP_Constraint(12, 1, 0.0, true);
	SP_Constraint *sp6 = new SP_Constraint(12, 2, 0.0, true);
	SP_Constraint *sp7 = new SP_Constraint(22, 0, 0.0, true);
	SP_Constraint *sp8 = new SP_Constraint(22, 1, 0.0, true);
	SP_Constraint *sp9 = new SP_Constraint(22, 2, 0.0, true);

	theDomain->addSP_Constraint(sp4);
	theDomain->addSP_Constraint(sp5);
	theDomain->addSP_Constraint(sp6);
	theDomain->addSP_Constraint(sp7);
	theDomain->addSP_Constraint(sp8);
	theDomain->addSP_Constraint(sp9);

	Truss *truss1 = new Truss(12, 2, 8, 12, *elastic, 1.0);
	Truss *truss2 = new Truss(22, 2, 8, 22, *elastic, 1.0);

	theDomain->addElement(truss1);
	theDomain->addElement(truss2);

	//////////// add vertical loads

	TimeSeries *theSeries1 = new LinearSeries();

	LoadPattern *theLoadPattern1 = new LoadPattern(1);

	theLoadPattern1->setTimeSeries(theSeries1);
	theDomain->addLoadPattern(theLoadPattern1);

	Vector theLoadValues1(3);
	theLoadValues1(1) = 1.0e12;
	NodalLoad *theLoad1 = new NodalLoad(1, 8, theLoadValues1);
	theDomain->addNodalLoad(theLoad1, 1);

	//////////////  gravity analysis

	double incrFacotr1 = -1.0e-3;
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

	theAnalysis1.analyze(200);
	analysisStep += 200;

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

	double incrFacotr2 = 1.0e-3;

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

	ofstream file1, file2, file3, file4;
	file1.precision(12);
	file2.precision(12);
	file3.precision(12);
	file4.precision(12);

	/////////////////////////////////

	int numEle = 2;

	FiberSectionTimo **theSec = new FiberSectionTimo *[numEle];
	theSec[0] = (FiberSectionTimo *)beam1->theSection;
	//theSec[1] = (FiberSectionTimo *)beam2->theSection;
	//theSec[2] = (FiberSectionTimo *)beam3->theSection;
	theSec[1] = (FiberSectionTimo *)beam4->theSection;
	
	for (int count = 0; count < 2; count++)
	{
		std::string s1, s2, s3, s4;

		s1 = "D:\\pertSecShearForce" + std::to_string(count) + ".out";
		s2 = "D:\\pertSecMoment" + std::to_string(count) + ".out";
		s3 = "D:\\pertSecShearDeform" + std::to_string(count) + ".in";
		s4 = "D:\\pertSecAngle" + std::to_string(count) + ".in";

		file1.open(s1);
		file2.open(s2);
		file3.open(s3);
		file4.open(s4);

		double pertDispValue = 1.0e-10;
		double pertRotateValue = 1.0e-12;

		Matrix perSecForce(2 * halfPertStep + 1, 2 * halfPertStep + 1);
		Matrix perSecMoment(2 * halfPertStep + 1, 2 * halfPertStep + 1);

		Vector secStrain = theSec[count]->getSectionDeformation();

		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			Vector secStrainTr(3), secForce;

			for (int j = 0; j < 2 * halfPertStep + 1; j++)
			{
				secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + j * pertRotateValue;
				secStrainTr(2) = secStrain(2) - halfPertStep * pertDispValue + i * pertDispValue;

				theSec[count]->setTrialSectionDeformation(secStrainTr);

				secForce = theSec[count]->getStressResultant();
				perSecForce(i, j) = secForce(2);
				perSecMoment(i, j) = secForce(1);

				file1 << perSecForce(i, j) << "\t";
				file2 << perSecMoment(i, j) << "\t";
				if (i == 0)
					file4 << secStrainTr(1) << endln;
			}
			file1 << endln;
			file2 << endln;
			file3 << secStrainTr(2) << endln;
		}
		file1.close();
		file2.close();
		file3.close();
		file4.close();
	}


	///////////////////////////////////



	//NDMaterial **matPtr = new NDMaterial *[16];
	//for (int i = 0; i < 4; i++)
	//	for (int j = 0; j < 4; j++)
	//		matPtr[i * 4 + j] = theSec[i]->theMaterials[j];

	//for (int count = 0; count < 16; count++)
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

	//	Matrix pertNormStress(2 * halfPertStep + 1, 2 * halfPertStep + 1);
	//	Matrix pertShearStress(2 * halfPertStep + 1, 2 * halfPertStep + 1);

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
	//			pertNormStress(i, j) = stress(0);
	//			pertShearStress(i, j) = stress(1);

	//			file1 << pertNormStress(i, j) << "\t";
	//			file2 << pertShearStress(i, j) << "\t";
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



	opserr << "finished." << endln;
	getchar();

	return 0;
}

