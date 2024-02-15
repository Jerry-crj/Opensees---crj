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

#ifndef ShearWall_h
#define ShearWall_h

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

class ShearWall : public Element {

public:
	ShearWall();
	ShearWall(int tag, int node1, int node2, int node3, int node4,
		SectionForceDeformation *theSections, BeamIntegration *bi);
	virtual ~ShearWall();

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

//	static const double one_over_root3;
//	static const double sg[2];
//	static const double wg[8];

	//static data
	static Matrix stiff;
	static Vector resid;
	static Matrix mass;
	static Matrix damping;

	//node information
	ID connectedExternalNodes;  //four node numbers
	Node *nodePointers[4];      //pointers to four nodes

	//material information
	enum { maxNumSections = 20 };
//	static const int numSectionsN = 1;	// number of fiber in cross section
//	static const int numSectionsS = 1;	// number of fiber in cross section
	static const int numSections = 1;	// number of fiber in cross section
	SectionForceDeformation *theSections;  //material pointer
	//SectionForceDeformation *theSectionsSy;  //material pointer
	//SectionForceDeformation *theSectionsSz;  //material pointer

	//local nodal coordinates, six coordinates for each of four nodes
//	static double xl[2][4];
	double H;
	double L;
	// vector for applying loads
	
	double b[3];		// Body forces
	double appliedB[6];		// Body forces applied with load
	int applyLoad;

	BeamIntegration *beamInt;
	Vector *load;
	Matrix *Ki;
	bool doUpdateBasis;

	void formResistingForce(void);
	void formTangentStiff(void);
	void formMass(void);
	void formInertiaLoad(void);

	void computeBasis(void);

	void getShapeFnc(double gaussCoor, double shp[2][2]);
	const Matrix& computeB(int i, double shp[2][2]);

	Matrix *crdTransf;
	Matrix *crdTransfTran;
	//static double workArea[];

};



#endif
