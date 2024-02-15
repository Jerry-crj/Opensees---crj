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
	remove("D://secInfoOfFailedStep.txt");
	remove("D://fiberInfoOfFailedStep.txt");

	int analysisStep = 0;

	Domain *theDomain = new Domain();

	//////////// build model

	/////////////////  beam
	Node *node1 = new Node(1, 3, 0.0, 0.0);
	Node *node2 = new Node(2, 3, 0.0, 2000.0);

	theDomain->addNode(node1);
	theDomain->addNode(node2);

	SP_Constraint *sp1 = new SP_Constraint(1, 0, 0.0, true);
	SP_Constraint *sp2 = new SP_Constraint(1, 1, 0.0, true);
	SP_Constraint *sp3 = new SP_Constraint(1, 2, 0.0, true);
	theDomain->addSP_Constraint(sp1);
	theDomain->addSP_Constraint(sp2);
	theDomain->addSP_Constraint(sp3);

	NDMaterial *concrete = new ovalConcrete(1, -40.0, -2.0e-3, -2.0e-2, 1.0e-4, 1.0e-3, 0.9, 0.7);
	NDMaterial *steel = new steelTimo(2, 336.0, 2.0e5, 0.035, 18, 0.925, 0.15, 0.1);
	UniaxialMaterial *elastic = new ElasticMaterial(3, 1.0e18);

	Fiber **fibers = new Fiber*[7];
	fibers[0] = new NDFiber2d(1, *concrete, 30294.0, -396);
	fibers[1] = new NDFiber2d(2, *concrete, 30294.0, -198);
	fibers[2] = new NDFiber2d(3, *concrete, 30294.0, 0);
	fibers[3] = new NDFiber2d(4, *concrete, 30294.0, 198);
	fibers[4] = new NDFiber2d(5, *concrete, 30294.0, 396);
	fibers[5] = new NDFiber2d(6, *steel, 1300.0, 440.0);
	fibers[6] = new NDFiber2d(7, *steel, 1300.0, -440.0);

	SectionForceDeformation* theSection = new FiberSectionTimo(1, 7, fibers);

	TimoshenkoModified* beam1 = new TimoshenkoModified(1, 1, 2, theSection, 0.0);

	theDomain->addElement(beam1);

	/////////////////  truss
	Node *node101 = new Node(101, 3, 0.0, 2001.0);
	Node *node102 = new Node(102, 3, 1.0, 2000.0);

	theDomain->addNode(node101);
	theDomain->addNode(node102);

	SP_Constraint *sp4 = new SP_Constraint(101, 0, 0.0, true);
	SP_Constraint *sp5 = new SP_Constraint(101, 1, 0.0, true);
	SP_Constraint *sp6 = new SP_Constraint(101, 2, 0.0, true);
	SP_Constraint *sp7 = new SP_Constraint(102, 0, 0.0, true);
	SP_Constraint *sp8 = new SP_Constraint(102, 1, 0.0, true);
	SP_Constraint *sp9 = new SP_Constraint(102, 2, 0.0, true);

	theDomain->addSP_Constraint(sp4);
	theDomain->addSP_Constraint(sp5);
	theDomain->addSP_Constraint(sp6);
	theDomain->addSP_Constraint(sp7);
	theDomain->addSP_Constraint(sp8);
	theDomain->addSP_Constraint(sp9);

	Truss *truss1 = new Truss(101, 2, 2, 101, *elastic, 1.0);
	Truss *truss2 = new Truss(102, 2, 2, 102, *elastic, 1.0);

	theDomain->addElement(truss1);
	theDomain->addElement(truss2);

	//////////// add vertical loads

	TimeSeries *theSeries1 = new LinearSeries();

	LoadPattern *theLoadPattern1 = new LoadPattern(1);

	theLoadPattern1->setTimeSeries(theSeries1);
	theDomain->addLoadPattern(theLoadPattern1);

	Vector theLoadValues1(3);
	theLoadValues1(1) = 1.0e18;
	NodalLoad *theLoad1 = new NodalLoad(1, 2, theLoadValues1);
	theDomain->addNodalLoad(theLoad1, 1);

	////////////  gravity analysis

	double incrFacotr1 = -1.0e-5;
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

	/// const load

	theDomain->setLoadConstant();
	theDomain->setCurrentTime(0.0);
	theDomain->setCommittedTime(0.0);

	//////////// add push over loads

	TimeSeries *theSeries2 = new LinearSeries();

	LoadPattern *theLoadPattern2 = new LoadPattern(2);

	theLoadPattern2->setTimeSeries(theSeries2);
	theDomain->addLoadPattern(theLoadPattern2);

	Vector theLoadValues2(3);
	theLoadValues2(0) = 1.0e18;
	NodalLoad *theLoad2 = new NodalLoad(2, 2, theLoadValues2);

	theDomain->addNodalLoad(theLoad2, 2);

	////////////  push over analysis

	double incrFacotr2 = 1.0e-4;

	int totalIterNum = 8;
	Vector maxIterDisp(totalIterNum);
	maxIterDisp(0) = 3.788;
	maxIterDisp(1) = 9.986;
	maxIterDisp(2) = 16.41;
	maxIterDisp(3) = 24.22;
	maxIterDisp(4) = 38.91;
	maxIterDisp(5) = 54.52;
	maxIterDisp(6) = 71.05;
	maxIterDisp(7) = 71.28;

	ofstream file;
	file.open("info.txt");

	bool isContinue =  true;
	for (int i = 0; i < totalIterNum; i++)
	{
		if (isContinue)
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
				{
					isContinue = false;
					break;
				}
				count++;
				analysisStep++;

				opserr << analysisStep << endln;
			}
		}
		else 
		{
			break;
		}
	}
	file.close();

	///// perturb

	int halfPertStep = 100;

	int numEle = 1;
	int numFibers = 5;

	double pertDispValue = 1.0e-7;
	double pertRotateValue = 1.0e-9;

	//TimoshenkoModified **theEle = new TimoshenkoModified *[numEle];
	//theEle[0] = beam1;

	FiberSectionTimo **theSec = new FiberSectionTimo *[numEle];
	theSec[0] = (FiberSectionTimo *)beam1->theSection;

	ofstream file1, file2, file3;

	file1.precision(12);
	file2.precision(12);
	file3.precision(12);

	for (int count = 0; count < numEle; count++)
		theSec[count]->revertToLastCommit();

	for (int count = 0; count < numEle; count++)
	{

		std::string s1;

		s1 = "D:\\pertDisp" + std::to_string(count + 1) + ".in";

		file1.open(s1);

		Vector secStrain = theSec[count]->getSectionDeformation();

		Vector secStrainTr(3), secForce;
		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			secStrainTr(0) = secStrain(0) - halfPertStep * pertDispValue + i * pertDispValue;
			secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + i * pertRotateValue;
			secStrainTr(2) = secStrain(2) - halfPertStep * pertDispValue + i * pertDispValue;

			file1 << secStrainTr(0) << "\t" << secStrainTr(1) << "\t" << secStrainTr(2) << "\t" << endln;
		}
		file1.close();
	}

	for (int count = 0; count < numEle; count++)
		theSec[count]->revertToLastCommit();

	for (int count = 0; count < numEle; count++)
	{
		std::string s1, s2, s3;

		s1 = "D:\\normFPertDvDh" + std::to_string(count + 1) + ".out";
		s2 = "D:\\momentPertDvDh" + std::to_string(count + 1) + ".out";
		s3 = "D:\\shearFPertDvDh" + std::to_string(count + 1) + ".out";

		file1.open(s1);
		file2.open(s2);
		file3.open(s3);

		Vector secStrain = theSec[count]->getSectionDeformation();

		Vector secStrainTr(3), secForce;

		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			for (int j = 0; j < 2 * halfPertStep + 1; j++)
			{
					secStrainTr(0) = secStrain(0) - halfPertStep * pertDispValue + i * pertDispValue;
					secStrainTr(1) = secStrain(1);
					secStrainTr(2) = secStrain(2) - halfPertStep * pertDispValue + j * pertDispValue;

					theSec[count]->setTrialSectionDeformation(secStrainTr);

					secForce = theSec[count]->getStressResultant();

					file1 << secForce(0) << "\t";
					file2 << secForce(1) << "\t";
					file3 << secForce(2) << "\t";
			}
			file1 << endln;
			file2 << endln;
			file3 << endln;

		}

		file1.close();
		file2.close();
		file3.close();
	}

	for (int count = 0; count < numEle; count++)
		theSec[count]->revertToLastCommit();

	for (int count = 0; count < numEle; count++)
	{
		std::string s1, s2, s3;

		s1 = "D:\\normFPertDvDt" + std::to_string(count + 1) + ".out";
		s2 = "D:\\momentPertDvDt" + std::to_string(count + 1) + ".out";
		s3 = "D:\\shearFPertDvDt" + std::to_string(count + 1) + ".out";

		file1.open(s1);
		file2.open(s2);
		file3.open(s3);

		Vector secStrain = theSec[count]->getSectionDeformation();

		Vector secStrainTr(3), secForce;

		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			for (int j = 0; j < 2 * halfPertStep + 1; j++)
			{
				secStrainTr(0) = secStrain(0) - halfPertStep * pertDispValue + i * pertDispValue;
				secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + j * pertRotateValue;
				secStrainTr(2) = secStrain(2);

				theSec[count]->setTrialSectionDeformation(secStrainTr);

				secForce = theSec[count]->getStressResultant();

				file1 << secForce(0) << "\t";
				file2 << secForce(1) << "\t";
				file3 << secForce(2) << "\t";
			}
			file1 << endln;
			file2 << endln;
			file3 << endln;

		}
		file1.close();
		file2.close();
		file3.close();
	}

	for (int count = 0; count < numEle; count++)
		theSec[count]->revertToLastCommit();

	for (int count = 0; count < numEle; count++)
	{
		std::string s1, s2, s3;

		s1 = "D:\\normFPertDhDt" + std::to_string(count + 1) + ".out";
		s2 = "D:\\momentPertDhDt" + std::to_string(count + 1) + ".out";
		s3 = "D:\\shearFPertDhDt" + std::to_string(count + 1) + ".out";

		file1.open(s1);
		file2.open(s2);
		file3.open(s3);

		Vector secStrain = theSec[count]->getSectionDeformation();

		Vector secStrainTr(3), secForce;

		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			for (int j = 0; j < 2 * halfPertStep + 1; j++)
			{
				secStrainTr(0) = secStrain(0);
				secStrainTr(1) = secStrain(1) - halfPertStep * pertRotateValue + i * pertRotateValue;
				secStrainTr(2) = secStrain(2) - halfPertStep * pertDispValue + j * pertDispValue;

				theSec[count]->setTrialSectionDeformation(secStrainTr);

				secForce = theSec[count]->getStressResultant();

				file1 << secForce(0) << "\t";
				file2 << secForce(1) << "\t";
				file3 << secForce(2) << "\t";
			}
			file1 << endln;
			file2 << endln;
			file3 << endln;

		}

		file1.close();
		file2.close();
		file3.close();
	}


	///////////////////////////////////

	double pertNormStrain = 2.e-7;
	double pertShearStrain = 2.e-7;

	for (int count = 0; count < numEle; count++)
		theSec[count]->revertToLastCommit();

	NDMaterial **matPtr = new NDMaterial *[5];
	for (int j = 0; j < numFibers; j++)
		matPtr[j] = theSec[0]->theMaterials[j];

	for (int count = 0; count < numEle * numFibers; count++)
	{
		std::string s1;

		s1 = "D:\\pertStrain" + std::to_string(count + 1) + ".in";

		file1.open(s1);

		Vector strain = matPtr[count]->getStrain();

		Vector strainTr(2);

		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
			strainTr(0) = strain(0) - halfPertStep * pertNormStrain + i * pertNormStrain;
			strainTr(1) = strain(1) - halfPertStep * pertShearStrain + i * pertShearStrain;
			file1 << strainTr(0) << "\t" << strainTr(1) << endln;
		}
		file1.close();
	}

	for (int count = 0; count < numEle * numFibers; count++)
	{
		std::string s1, s2, s3, s4;

		s1 = "D:\\normStressPert" + std::to_string(count + 1) + ".out";
		s2 = "D:\\shearStressPert" + std::to_string(count + 1) + ".out";

		file1.open(s1);
		file2.open(s2);

		Vector strain = matPtr[count]->getStrain();

		Vector strainTr(2), stress(2);

		for (int i = 0; i < 2 * halfPertStep + 1; i++)
		{
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

