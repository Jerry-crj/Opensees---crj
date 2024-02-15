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
                                                                        
// $Revision: 1.3 $
// $Date: 2006-11-03 18:40:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/CombinedUniaxialTimoMat.h,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class definition for 
// CombinedUniaxialTimoMat. CombinedUniaxialTimoMat is based on an f2c of the FEDEAS material
// CombinedUniaxialTimoMat.f which is:
//-----------------------------------------------------------------------
// MENEGOTTO-PINTO STEEL MODEL WITH FILIPPOU ISOTROPIC HARDENING
//            written by MOHD YASSIN (1993)
//          adapted to FEDEAS material library
//    by E. Spacone, G. Monti and F.C. Filippou (1994)
//-----------------------------------------------------------------------


#ifndef CombinedUniaxialTimoMat_h
#define CombinedUniaxialTimoMat_h

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class CombinedUniaxialTimoMat : public NDMaterial
{
  public:
    CombinedUniaxialTimoMat(int tag, UniaxialMaterial* mat1, UniaxialMaterial* mat2);
    	    
    CombinedUniaxialTimoMat(void);
    virtual ~CombinedUniaxialTimoMat();
    

    const char *getType(void) const {return "CombinedUniaxialTimoMat";};

    NDMaterial *getCopy(void);
	NDMaterial *getCopy(const char *type);

    int setTrialStrain(const Vector &strainFromElement);
	
	//send back the strain
	const Vector& getStrain();

	//send back the stress 
	const Vector& getStress();

	//send back the tangent 
	const Matrix& getTangent();

	const Matrix& getInitialTangent();

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int	getResponse(int responseID, Information& matInfo);

 protected:
    
 private:
     UniaxialMaterial** getMaterials();

     UniaxialMaterial** theMaterials;
     Vector sig, sigP, eps, epsP;
     Matrix k, kP;

     Matrix* ki;
};


#endif

