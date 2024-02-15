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
                                                                        
// $Revision: 1.14 $
// $Date: 2008-08-26 16:47:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionTimoQualShear.h,v $
                                                                        
// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for 
// FiberSectionTimoQualShear.h. FiberSectionTimoQualShear provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSectionTimoQualShear_h
#define FiberSectionTimoQualShear_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

#define SEC_TAG_FiberSectionTimoQualShear 8532151

class NDMaterial;
class Fiber;
class Response;
class SectionIntegration;

class FiberSectionTimoQualShear : public SectionForceDeformation
{
  public:
    FiberSectionTimoQualShear(); 
    FiberSectionTimoQualShear(int tag, int numFibers, Fiber **fibs );
    ~FiberSectionTimoQualShear();

    const char *getClassType(void) const {return "FiberSectionTimoQualShear";};

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
 
    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int addFiber(Fiber &theFiber);

    void setSectionSize(double b0, double h0);

  protected:
    
  private:

    int numFibers, sizeFibers;       // number of fibers in the section
    NDMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc, zloc, area]

    double QzBar, Abar;
    double yBar;      // Section centroid
    double b, h;

    static ID code;
    SectionIntegration *sectionIntegr;

    Vector e;          // trial section deformations 
	Vector s;         // section resisting forces  (axial force, bending moment)
    Matrix ks;        // section stiffness
	Matrix *ki;

};

#endif
