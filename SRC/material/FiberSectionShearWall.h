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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionShearWall.h,v $
                                                                        
// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for 
// FiberSectionShearWall.h. FiberSectionShearWall provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSectionShearWall_h
#define FiberSectionShearWall_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

#define SEC_TAG_FiberSectionShearWall 100000001

class NDMaterial;
class Fiber;
class Response;
class SectionIntegration;

class FiberSectionShearWall : public SectionForceDeformation
{
  public:
    FiberSectionShearWall(); 
    FiberSectionShearWall(int tag, int numFibers, Fiber **fibs );
    ~FiberSectionShearWall();

    const char *getClassType(void) const {return "FiberSectionShearWall";};

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

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);

    //const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    //const Matrix & getSectionTangentSensitivity(int gradIndex);
    //int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    //const Vector & getSectionDeformationSensitivity(int gradIndex);
    // AddingSensitivity:END ///////////////////////////////////////////

	//by SAJalali
	//double getEnergy() const;


  protected:
    
  private:
    int numFibers, sizeFibers;       // number of fibers in the section
    NDMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc, zloc, area]
    double   kData[64];              // data for ks matrix 
    double   sData[8];               // data for s vector 

    double QzBar, QyBar, Abar;
    double yBar;       // Section centroid
    double zBar;
  
    SectionIntegration *sectionIntegr;

    static ID code;

    Vector e;          // trial section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness

};

#endif
