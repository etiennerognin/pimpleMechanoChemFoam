/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    pimpleMechanoChemFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    
    #include "createMechanoChemFields.H"
    
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    //turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    fileName outputFile(runTime.path()/"time_conversion_pressure.txt");
    OFstream os(outputFile);
    
    // Momentum predictor at first iteration
    bool firstIter = true;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
            
            
            #include "AEqn.H"
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            
            
            
            
            /*
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
            */
            if (not firstIter)
            {
                #include "cEqn.H"
            }
            
            
        }
        firstIter = false;

        
        if (runTime.outputTime())
        {
            // Record file

            // Calculate output conversion flux:
            
            // Flux at outlet
            label patchID(mesh.boundaryMesh().findPatchID("outlet"));
            dimensionedScalar convFlux1 = gSum(conv.boundaryField()[patchID]*phi.boundaryField()[patchID])/gSum(phi.boundaryField()[patchID]);
            
            // Volume integral of the reaction rate
            dimensionedScalar convFlux2 = fvc::domainIntegrate(kc)/gSum(-phi.boundaryField()[mesh.boundaryMesh().findPatchID("inlet")]);
            
            // Average the input pressure
            patchID = mesh.boundaryMesh().findPatchID("inlet");
            dimensionedScalar reducedP = gSum(p.boundaryField()[patchID]*mesh.magSf().boundaryField()[patchID])/gSum(mesh.magSf().boundaryField()[patchID])*rho*1e-5;
            
            os << runTime.timeName() << "\t" << convFlux1.value() << "\t" << convFlux2.value() << "\t" << reducedP.value() << endl;
        }
        
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
