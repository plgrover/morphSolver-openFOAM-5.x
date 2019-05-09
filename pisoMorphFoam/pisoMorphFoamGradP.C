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
#include "turbulenceModel.H"
#include "RASModel.H"
// #include <math.h>       /* atan */
#define PI 3.14159265


#include "Avalanche.H"
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
    // #include "Avalanche.H"    

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
	// https://www.cfd-online.com/Forums/openfoam-programming-development/122557-moving-boundary-problem-based-calculated-data.html#post447552
	//Find the relevant size of the vector and declare a vectorField.
	int Psize= pDisp.size();
	vectorField dispVals(Psize);

	Info<< "PSize: " << Psize << nl << endl;

	vectorField &PointPointer = refCast<vectorField>(PointDisplacement.boundaryFieldRef()[patchWallID]);
	vectorField PointNormalVector = mesh.boundaryMesh()[patchWallID].pointNormals();
	
	// Get ahold of the turbulent properties
	const tmp<scalarField> tnuw = turbulence().nu(patchWallID);
	const scalarField& nuw = tnuw();
	// Info<< "Nu: " << nuw << nl << endl;
	
	const tmp<scalarField> tnutw = turbulence().nut(patchWallID);         
	const scalarField& nutw = tnutw();
	// Info<< "Nu: " << nutw << nl << endl;
	
	const fvPatchVectorField& Uw = turbulence().U().boundaryField()[patchWallID];
    	const scalarField magGradUw((mag(Uw.snGrad())));

	//const scalarField tau = magGradUw*(nutw + nuw);
	const scalar D50 = 0.001;
	const scalar gammaS = 9.81*(2650. - 1000.);
	
	
	scalarField Y = (1000.)*magGradUw*(nutw + nuw)/scalar(D50*gammaS);
	
	scalarField tau = (nuw + nutw)*magGradUw*scalar(1000.);

	//Info << "Tau " << tau << endl;

	const scalar Ycr = 0.6/(D50*gammaS);
			

	scalarField qbs = (0.5*8.5)*(Foam::sqrt(mag(Y)))*(Y - Ycr); 
	// Get access to the polyPatch

	// Info << "qbs: " << qbs << endl; 

	const polyPatch& pp = mesh.boundaryMesh()[patchWallID];
	
	scalarField dzPatch =  Y * (0);

	forAll(qbs, index)
	{
		double qbi = qbs[index]; 
		double qbim1 = 0.0;
		double xi = pp.faceCentres()[index].x();
		double xim1 = 0.0; 
		
		
		if (index > 0)
		{
			qbim1 = qbs[index - 1];
			xim1 = pp.faceCentres()[index].x() ;
			xi = pp.faceCentres()[index + 1].x();
		}
		else
		{
			qbim1 = qbs[0];
		}
		
		if (qbi < 0.) 
		{
			qbi = 0.;
		}
		if (qbim1 < 0.)
		{
			qbim1 = 0.;
		}

		double dx = xi - xim1;
		double dq = qbi - qbim1;
		if (index < 2)
		{
			dzPatch[index] = 0.;
		}
		else
		{
			dzPatch[index] = runTime.deltaT().value()/(2.*dx*(1.-0.4))*dq;
		}
		// Info << "Ycr " << Ycr << "  Y: " << Y[index] << " dq " << dq << " dz: " << dzPatch[index] << endl;		
	}

	Info << "dz size " << dzPatch.size() << endl;

	primitivePatchInterpolation patchInterpolator (mesh.boundaryMesh()[patchWallID]);
	scalarField facedZValues = patchInterpolator.faceToPointInterpolate(dzPatch);
	Info << "facedzvalues size " << facedZValues.size() << endl;

	// int nuSize = nuw.size();
	//scalarField<double> qbed(nuSize);

	Info << "qobs " << qbs.size() << endl;
	// pp size is equal to the number of faces. 
	Info<< "pp size: " << pp.size() << endl;
	
	// Gets the list of points on the patch 
	//const List<vector>& patchFound = mesh.boundaryMesh()[patchWallID].localPoints(); 

    // Testing Avalanche here
	vectorField test = pp.faceCentres();
	// Info << "Type of " << TypeNameNoDebug(pp.faceCentres()) << endl;
	avalancheProfile(test, 15.0, 15.2, 15.4); // pp.faceCentres, 30., 30., 30.);

 	forAll(dispVals, index)
	{	
		/*if (facedZValues[index] > 0.0)
		{
			Info<< "Dz: " << facedZValues[index] << endl;
		}*/
		dispVals[index].x() = PointPointer[index].x();
		dispVals[index].y() = PointPointer[index].y(); //- 100. * PointNormalVector[index].y() * runTime.deltaT().value();
		dispVals[index].z() = PointPointer[index].z() - facedZValues[index];
		// Info << "PP x: " << patchFound[index].x() << endl;
		 
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
