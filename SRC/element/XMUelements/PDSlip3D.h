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

// $Revision: 1.15 $
// $Date: 2009-08-07 20:01:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PDSlip/PDSlip.h,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for PDSlip.

#ifndef PDSlip3D_h
#define PDSlip3D_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#define ELE_TAG_PDSlip3D  123123111

class Node;
class NDMaterial;
class UniaxialMaterial;
class Response;

class PDSlip3D : public Element
{
public:
	PDSlip3D(int tag, int nd1, int nd2, int nd3, int nd4, int nd5,
		int nd6, int nd7, int nd8, int nd9,
		UniaxialMaterial& m1, UniaxialMaterial& m2,
		double p1, double p2, double p3);
	PDSlip3D();
	~PDSlip3D();

	const char* getClassType(void) const { return "PDSlip3D"; }

	int getNumExternalNodes(void) const;
	const ID& getExternalNodes(void);
	Node** getNodePtrs(void);

	int getNumDOF(void);
	void setDomain(Domain* theDomain);

	// public methods to set the state of the element    
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	int update(void);

	// public methods to obtain stiffness, mass, damping and residual information    
	const Matrix& getTangentStiff(void);
	const Matrix& getInitialStiff(void);
	const Matrix& getMass(void);

	void zeroLoad();
	int addLoad(ElementalLoad* theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector& accel);

	const Vector& getResistingForce(void);
	const Vector& getResistingForceIncInertia(void);

	// public methods for element output
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker
		& theBroker);

	int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
	void Print(OPS_Stream& s, int flag = 0);

	Response* setResponse(const char** argv, int argc,
		OPS_Stream& s);

	int getResponse(int responseID, Information& eleInformation);

private:

	void determinPst();

	void shapeFunction(double xi, double eta, double zeta);

	void getTransf();

	//const Vector& getResistingForce(Vector eps);

	// private attributes - a copy for each object of the class

	UniaxialMaterial** theSprMaterial; // pointer to the ND material objects

	ID connectedExternalNodes; // Tags of quad nodes

	Node* theNodes[9];

	static double matrixData[729];  // array data for matrix
	static Matrix K;		// Element stiffness, damping, and mass Matrix
	static Vector P;		// Element resisting force vector

	static double shp[8];	// Stores shape functions and derivatives (overwritten)

	// private member functions - only objects of this class can call these

	Matrix* Ki;
	Vector dv;

	Vector pst;
	Matrix T, Ttran;

	double deform, force;
};

#endif

