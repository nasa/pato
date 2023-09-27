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

Class
    Foam::simpleFoam

Description

SourceFiles
    simpleFoam.H
    simpleFoam.C

\*---------------------------------------------------------------------------*/

#ifndef simpleFoam_H
#define simpleFoam_H

#include "basicFluidModel.H"
// ADDED INCLUDE
#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class simpleFoam
:
public basicFluidModel
{

 protected:
  // createFluidMeshes
  PtrList<dynamicFvMesh> fluidRegions;
  // createFluidFields
  PtrList<volVectorField> UFluid;
  PtrList<surfaceScalarField> phiFluid;
  PtrList<volScalarField> pFluid;
  PtrList<incompressible::turbulenceModel> turbulence;
  PtrList<singlePhaseTransportModel> transport;
  PtrList<IOMRFZoneList> MRFfluid;
  PtrList<fv::options> fluidFvOptions;
  scalarList cumulativeContErr;

 public:

  //- Runtime type information
  TypeName("simpleFoam");

  // Constructors

  //- Construct from mesh and phase name
  simpleFoam(Time& runTime);

  //- Destructord
  virtual ~simpleFoam();

  //- Member Functions

  //- Update before runTime++ in the while loop
  virtual void updateBefore();

  //- Update after runTime++ in the while loop
  virtual void updateAfter();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "simpleFoam.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
