/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoamGradP

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "pointMesh.H"
#include "dynamicFvMesh.H"
#include "primitivePatchInterpolation.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    //#include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
	// - Update of the mesh goes here
	// Get the patch ID where impingement occurs: patch name must be "wall"
    	label patchWallID = mesh.boundaryMesh().findPatchID("bottomWall");
    	const fvPatch& patchWallFaces = mesh.boundary()[patchWallID];

    	Info<< "#--------------- Displacement ---------------# " << endl;

	//Find the reference to the location of pointDisplacement field
	pointVectorField& PointDisplacement = const_cast<pointVectorField&>
	(
		mesh.objectRegistry::lookupObject<pointVectorField>
		(
			"pointDisplacement"
		)
	);
	
	//Get the vector field of the patch
	vectorField &pDisp=refCast<vectorField>(PointDisplacement.boundaryFieldRef()[patchWallID]);

	//Find the relevant size of the vector and declare a vectorField.
	int Psize= pDisp.size();
	vectorField dispVals(Psize);

	Info<< "PSize: " << Psize << nl << endl;

	//- set-up interpolator
	primitivePatchInterpolation patchInterpolator (mesh.boundaryMesh()[patchWallID] );

	/*scalarField MimpPatch = Mimp.boundaryField()[patchWallID];

	//- perform interpolation
	scalarField faceValues = patchInterpolator.faceToPointInterpolate(MimpPatch); */

	vectorField &PointPointer = refCast<vectorField>(PointDisplacement.boundaryFieldRef()[patchWallID]);
	vectorField PointNormalVector = mesh.boundaryMesh()[patchWallID].pointNormals();

 	forAll(dispVals, index)
	{
		dispVals[index].x() = PointPointer[index].x();
		dispVals[index].y() = PointPointer[index].y(); //- 100. * PointNormalVector[index].y() * runTime.deltaT().value();
		dispVals[index].z() = PointPointer[index].z() - 0.01 * runTime.deltaT().value();
	} 
	
	PointDisplacement.boundaryFieldRef()[patchWallID] == dispVals;	

	
	// - End of update meths 
	mesh.update();
        
	Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
