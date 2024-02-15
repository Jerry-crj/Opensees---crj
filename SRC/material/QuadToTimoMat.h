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
// $Date: 2007-06-08 00:38:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/QuadToTimoMat.h,v $
                                                                      
// Written: fmk
// Created: 03/06
//
// Description: This file contains the class definition for 
// QuadToTimoMat. QuadToTimoMat is based on an f2c of the FEDEAS material
// Concr2.f which is:
/*-----------------------------------------------------------------------
! concrete model with damage modulus    
!       by MOHD YASSIN (1993)
! adapted to FEDEAS material library
! by D. Sze and Filip C. Filippou in 1994
-----------------------------------------------------------------------*/



#ifndef QuadToTimoMat_h
#define QuadToTimoMat_h

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicThreeDimensional.h>

class QuadToTimoMat : public NDMaterial
{
  public:
    QuadToTimoMat(int tag, NDMaterial * matTag);

    QuadToTimoMat(void);

    ~QuadToTimoMat();

    const char *getType(void) const {return "QuadToTimoMat";};
    NDMaterial *getCopy(void);
	NDMaterial *getCopy(const char *type);
 
	int setTrialStrain(const Vector &v);
	int setTrialStrain(const Vector &v, const Vector &r);
	int setTrialStrainIncr(const Vector &v);
	int setTrialStrainIncr(const Vector &v, const Vector &r);
	const Matrix &getTangent(void);
	const Matrix &getInitialTangent(void) { return this->getTangent(); };

	const Vector &getStress(void);
	const Vector &getStrain(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
 
	Response * setResponse(const char **argv, int argc, OPS_Stream &output);
	int	getResponse(int responseID, Information &matInfo);

 protected:
    
 private:
    //void Tens_Envlp (double epsc, double &sigc, double &Ect);
    //void Compr_Envlp (double epsc, double &sigc, double &Ect);
	
	// ElasticIsotropicPlaneStrain2D * theMat;
	 NDMaterial * theMat;

	// stiffness without damaged occured
	 Vector stress;
	 Vector strain;
	 Matrix tangent;
};


#endif

