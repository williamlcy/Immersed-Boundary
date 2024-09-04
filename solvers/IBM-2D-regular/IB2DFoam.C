/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/
#include "mpi.h"
#include "fvCFD.H"
#include "pisoControl.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "baseFunction.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //make sure mpi environment is initialised correctly
int mpi_initised=-1;
MPI_Initialized(&mpi_initised);
//  调用 MPI_Initialized 函数检查 MPI 是否已经初始化。
//  如果已初始化，则 mpi_initised 会被设置为非零值；
//  如果未初始化，则 mpi_initised 保持为零
if(mpi_initised)
{
}
else
{
MPI_Init(&argc,&argv);
//  如果 MPI 未初始化（mpi_initised 为零），则调用 MPI_Init(&argc, &argv) 来初始化 MPI 环境。
}

    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "createIBFields.H"
    #include "initContinuityErrs.H"
    


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        #include "inputData.H"
        std::ifstream para_file("Paralist.txt");
        std::vector<vector> para_list;
        scalar para_X, para_Y, para_Z;
        if (para_file.is_open()) {

            Info<<"reading the file paralists.txt" << endl;
            while (para_file >> para_X >> para_Y >> para_Z) {
                vector DataTransfer(para_X, para_Y, para_Z);
                para_list.push_back(DataTransfer);
            }

            file.close();
        } else{
            Info<<"can not find the file" << endl;
            para_X = 0.0;
            para_Y = 0.0;
            para_Z = 0.0;
            for(int i = 0;i<runTime.endTime().value()/runTime.deltaT().value();i++)
            {
                vector DataTransfer(para_X, para_Y, para_Z);
                para_list.push_back(DataTransfer);
            }
        }
        if(para_list.size() != runTime.endTime().value()/runTime.deltaT().value())
        {
            FatalError << "The number of parameters is not equal to the number of time steps!" << endl;
        } else{
            label list_x = runTime.value()/runTime.deltaT().value()-1;
            Info << list_x << endl;
            for(label I=0;I<Lmarks.size();I++)
            {
                Lmarks[I].x() += para_list[list_x].x();
                Lmarks[I].y() += para_list[list_x].y();
                Lmarks[I].z() += para_list[list_x].z();
                Info << Lmarks[I] << endl;
            }

        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "meshSupportSearch.H"
        #include "CourantNo.H"

        // Momentum predictor (in this stage, the predicted U field does not satisfy the continuity equation)
    
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));  
        }

        #include "Interpolation.H"
        #include "Spreading.H"

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA) + rAU*fvc::div(F_Fluid)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA + dt*(linearInterpolate(F_Fluid) & mesh.Sf()) - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p) + rAU*F_Fluid;
            U.correctBoundaryConditions();
        }
        
        //scalar ET = U.mesh().time().value();

        #include "liftCoeff2.H"
        
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
