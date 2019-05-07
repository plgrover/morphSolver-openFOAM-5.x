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
	// https://www.cfd-online.com/Forums/openfoam-programming-development/122557-moving-boundary-problem-based-calculated-data.html#post447552
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
	
	// Get ahold of the turbulent properties
	// Extracted from wallShearStress.C code (v2.3)
	/*autoPtr<incompressible::RASModel> model     
	(         
		incompressible::RASModel::New(U, phi, laminarTransport)     
	);*/

	// const volScalarField& k = lookupObject<volScalarField>("k");
	const volScalarField& k = turbulence().k();
	const tmp<scalarField> tnuw = turbulence().nu(patchWallID);
	const scalarField& nuw = tnuw();
	Info<< "Nu: " << nuw << nl << endl;
	
	const tmp<scalarField> tnutw = turbulence().nut(patchWallID);         
	const scalarField& nutw = tnutw();
	Info<< "Nu: " << nutw << nl << endl;
	
	const fvPatchVectorField& Uw = turbulence().U().boundaryField()[patchWallID];
    	const scalarField magGradUw((mag(Uw.snGrad())));

	const scalarField tau = magGradUw*(nutw + nuw);
	scalarField qbs = Foam::sqrt(Foam::cmptMag(tau))*(tau - scalar(0.6));
	
	// Get access to the polyPatch
	const polyPatch& pp = mesh.boundaryMesh()[patchWallID];

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
		double dz = runTime.deltaT().value()/(2.*dx*(1.-0.4))*dq;
		
	}

	// int nuSize = nuw.size();
	//scalarField<double> qbed(nuSize);

	Info << "qobs " << qbs.size() << endl;
	// pp size is equal to the number of faces. 
	Info<< "pp size: " << pp.size() << endl;
	
	// Gets the list of points on the patch 
	const List<vector>& patchFound = mesh.boundaryMesh()[patchWallID].localPoints();

 	forAll(dispVals, index)
	{
		//vector pointX(0,0,0);
		// vector pos( mesh.Cf().boundaryField()[patchWallID][index] );
		double x = pp.faceCentres()[index].x();
		
		
		double qbUpwind = 0.0;
		double qb = 0.0;
		double xUpwind = 0.0;
		
		if (qbs[index] > 0.0)
		{
			qb = qbs[index];
		}

		if (index > 0)
		{	
			qbUpwind = qbs[index - 1];
			if (qbUpwind < 0.)
			{
				qbUpwind = 0.;
			}
			// vector pos( mesh.Cf().boundaryField()[patchWallID][index] );
			xUpwind = 0.0;  //pos[0];
		}

		double dx = x - xUpwind;
		// Info<<"indx: " << index << " x "  << pp.faceCentres()[index+5].x() << " y " << pp.faceCentres()[index].y() << " z " << pp.faceCentres()[index].z() << endl;
		// Info<<"------------"<< nl << endl;
		double dq = qbUpwind - qb;
		
		double dz = 0.0; //(runTime.deltaT().value()/(2.*dx*(1-0.4)))*dq;
		

		dispVals[index].x() = PointPointer[index].x();
		dispVals[index].y() = PointPointer[index].y(); //- 100. * PointNormalVector[index].y() * runTime.deltaT().value();
		dispVals[index].z() = PointPointer[index].z() - 0.001;
		Info << "PP x: " << patchFound[index].x() << endl;
		 
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
