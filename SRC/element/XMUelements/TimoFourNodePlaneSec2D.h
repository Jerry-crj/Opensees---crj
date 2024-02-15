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

// $Revision: 1.19 $
// $Date: 2007-10-13 00:53:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/TimoFourNodePlaneSec2D.h,v $

// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for TimoFourNodePlaneSec2D.
// The element displacement field gives rise to constant axial strain, constant shearing and
// linear curvature.

#ifndef TimoFourNodePlaneSec2D_H
#define TimoFourNodePlaneSec2D_H

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BeamIntegration.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class Response;

#define ELE_TAG_TimoFourNodePlaneSec2D  95561812123334

class TimoFourNodePlaneSec2D : public Element
{
public:
	TimoFourNodePlaneSec2D(int tag, int nd1, int nd2, int nd3, int nd4,
		int numSections, SectionForceDeformation** s,
		BeamIntegration& bi);
	TimoFourNodePlaneSec2D();
	~TimoFourNodePlaneSec2D();

	const char* getClassType(void) const { return "TimoFourNodePlaneSec2D"; };

	int getNumExternalNodes(void) const;
	const ID& getExternalNodes(void);
	Node** getNodePtrs(void);

	int getNumDOF(void);
	void setDomain(Domain* theDomain);

	// public methods to set the state of the element    
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// public methods to obtain stiffness, mass, damping and residual information    
	int update(void);
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
	int displaySelf(Renderer& theViewer, int displayMode, float fact);
	void Print(OPS_Stream& s, int flag = 0);

	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& eleInfo);

protected:

private:
	void calculateB(double xi);

	int numSections;
	SectionForceDeformation** theSections; // pointer to the ND material objects

	BeamIntegration* beamInt;

	ID connectedExternalNodes; // Tags of quad nodes

	Node* theNodes[4];

	static Matrix K;		// Element stiffness, damping, and mass Matrix
	static Vector P;		// Element resisting force vector

	Matrix Ki;
	static Matrix B, Btran;

	enum { maxNumSections = 20 };

//	static double workArea[];

	double L;
	double sintheta;
	double costheta;

	Matrix T, Ttran;
};

#endif

