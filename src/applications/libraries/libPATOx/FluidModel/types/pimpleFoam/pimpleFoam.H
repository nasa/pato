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
    Foam::pimpleFoam

Description

SourceFiles
    pimpleFoam.H
    pimpleFoam.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleFoam_H
#define pimpleFoam_H

#include "basicFluidModel.H"
// ADDED INCLUDE
#include "turbulentTransportModel.H"
#include "courantNo.H"
#include "singlePhaseTransportModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class pimpleFoam
:
public basicFluidModel
{

 protected:
  // createFluidMeshes
  PtrList<dynamicFvMesh> fluidRegions;
  // createFluidFields
  PtrList<volVectorField> UFluid;
  List<autoPtr<surfaceVectorField>> Uf_auto_list;
  PtrList<surfaceScalarField> phiFluid;
  PtrList<volScalarField> pFluid;
  PtrList<incompressible::turbulenceModel> turbulence;
  PtrList<singlePhaseTransportModel> transport;
  PtrList<IOMRFZoneList> MRFfluid;
  PtrList<fv::options> fluidFvOptions;
  // PIMPLE controls
  label nOuterCorr;
  bool correctPhi;
  // Time Controls
  scalarList cumulativeContErr;
  // incompressible Courant Number
  scalar CoNum;
  scalar MeanCoNum;
  // createTimeControls
  bool adjustTimeStep;
  scalar maxCo;
  scalar maxDeltaT;
  // Mesh motion controls
  bool checkMeshCourantNo;
  label moveMeshOuterCorrectors;
  // Reference pressure
  label pRefCell;
  scalar pRefValue;

 public:

  //- Runtime type information
  TypeName("pimpleFoam");

  // Constructors

  //- Construct from mesh and phase name
  pimpleFoam(Time& runTime);

  //- Destructord
  virtual ~pimpleFoam();

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
#include "pimpleFoam.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
