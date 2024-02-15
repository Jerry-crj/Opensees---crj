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

// $Revision: 1.10 $
// $Date: 2011/03/10 22:51:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/ShellMITC4.h,v $

// Written: Leopoldo Tesser, Diego A. Talledo, Véronique Le Corvec
//
// Bathe MITC 4 four node shell element with membrane and drill
// Ref: Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984

#ifndef PDTrussWithShear2D_h
#define PDTrussWithShear2D_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <R3vectors.h>
#include <CrdTransf.h>
#include <BeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <NDMaterial.h>

#define ELE_TAG_PDTrussWithShear2D 7854688349

class PDTrussWithShear2D: public Element {

public:
	PDTrussWithShear2D();
	PDTrussWithShear2D(int tag, int node1, int node2, double a, NDMaterial* theMat);
	virtual ~PDTrussWithShear2D();

	void	setDomain(Domain *theDomain);
	int		getNumExternalNodes() const;
	const	ID &getExternalNodes();
	Node	**getNodePtrs();
	int		getNumDOF();
	int		commitState();
	int		revertToLastCommit();
	int		revertToStart();

	//print out element data
	void Print(OPS_Stream &s, int flag);

	//return stiffness matrix 
	int update(void);
	const Matrix &getTangentStiff();
	const Matrix &getInitialStiff();
	const Matrix &getMass();

	// methods for applying loads
	void zeroLoad(void);
	int addLoad(ElementalLoad *theLoad, double loadFactor);
	int addInertiaLoadToUnbalance(const Vector &accel);

	//get residual
	const Vector &getResistingForce();
	const Vector &getResistingForceIncInertia();

	// public methods for element output
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
		&theBroker);

	Response* setResponse(const char **argv, int argc, OPS_Stream &output);
	int getResponse(int responseID, Information &eleInfo);

	//plotting 
	int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);

private:

	//static data
	static Matrix stiff;
	static Vector resid;
	static Matrix mass;
	static Matrix damping;

	//node information
	ID connectedExternalNodes;  //four node numbers
	Node *nodePointers[2];      //pointers to four nodes

	NDMaterial *theMaterial;  //material pointer

	double A;
	double L;
	double initL;
	Vector initdVec,initnVec;

	Vector *load;
	Matrix *Ki;

	Matrix T, Ttran;
};



#endif
