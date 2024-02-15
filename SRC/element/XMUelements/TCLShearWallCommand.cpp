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

// $Revision: 1.8 $
// $Date: 2009-08-07 00:18:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ShearWall/TclShearWallCommand.cpp,v $

// Written: fmk 
// Created: 03/01
//
// What: "@(#) TclShearWallCommand.cpp, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <LegendreBeamIntegration.h>
#include "ShearWall.h"

#ifdef _FLShearWall
#include <FLShearWall.h>
#endif

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addShearWall(ClientData clientData, Tcl_Interp *interp, int argc,
	TCL_Char **argv, Domain*theTclDomain,
	TclModelBuilder *theTclBuilder, int eleArgStart)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed\n";
		return TCL_ERROR;
	}

	// check the number of arguments is correct
	if ((argc - eleArgStart) < 7) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element ShearWall eleTag? Node1? Node2? Node3? Node4? secTag? crdTag\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, Node1, Node2, Node3, Node4, numSec, secTag;

	if (Tcl_GetInt(interp, argv[1 + eleArgStart], &eleTag) != TCL_OK) {
		opserr << "WARNING invalid ShearWall eleTag" << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2 + eleArgStart], &Node1) != TCL_OK) {
		opserr << "WARNING invalid Node1\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[3 + eleArgStart], &Node2) != TCL_OK) {
		opserr << "WARNING invalid Node2\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[4 + eleArgStart], &Node3) != TCL_OK) {
		opserr << "WARNING invalid Node3\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[5 + eleArgStart], &Node4) != TCL_OK) {
		opserr << "WARNING invalid Node4\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[6 + eleArgStart], &numSec) != TCL_OK) {
		opserr << "WARNING invalid numSec\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (numSec <= 0) {
		opserr << "WARNING invalid numSec - cannot be <= 0 \n";
		opserr << argv[1] << " element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[7 + eleArgStart], &secTag) != TCL_OK) {
		opserr << "WARNING invalid matTag\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}
	
	SectionForceDeformation **sections = new SectionForceDeformation *[numSec];
	SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
	if (theSection == 0) {
		opserr << "WARNING section not found\n";
		opserr << "Section: " << secTag;
		opserr << argv[1] << " element: " << eleTag << endln;
		return TCL_ERROR;
	}
	for (int i = 0; i < numSec; i++)
		sections[i] = theSection;
	
	//if (Tcl_GetInt(interp, argv[8 + eleArgStart], &bi) != TCL_OK) {
	//	opserr << "WARNING invalid bi\n";
	//	opserr << "ShearWall element: " << eleTag << endln;
	//	return TCL_ERROR;
	//}
	// now create the ShearWall and add it to the Domain
	BeamIntegration *bi = new LegendreBeamIntegration();

	Element *theShearWall = 0;
	if (strcmp(argv[1], "ShearWall") == 0) {
		theShearWall = new ShearWall(eleTag, Node1, Node2, Node3, Node4,
			numSec, sections, bi);
	}

	delete[] sections;

	if (theShearWall == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "ShearWall element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theShearWall) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "ShearWall element: " << eleTag << endln;
		delete theShearWall;
		return TCL_ERROR;
	}

	// if get here we have successfully created the node and added it to the domain
	return TCL_OK;
}


