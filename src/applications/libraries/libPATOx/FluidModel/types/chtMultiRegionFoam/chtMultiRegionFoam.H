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
    Foam::chtMultiRegionFoam

Description

SourceFiles
    chtMultiRegionFoam.H
    chtMultiRegionFoam.C

\*---------------------------------------------------------------------------*/

#ifndef chtMultiRegionFoam_H
#define chtMultiRegionFoam_H

#include "basicFluidModel.H"
// ADDED INCLUDE
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "compressibleCourantNo.H"
#include "radiationModel.H"
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

class chtMultiRegionFoam
:
public basicFluidModel
{

 protected:
  // createFluidMeshes
  PtrList<dynamicFvMesh> fluidRegions;
  // createFluidFields
  PtrList<rhoThermo> thermoFluid;
  PtrList<volScalarField> rhoFluid;
  PtrList<volVectorField> UFluid;
  PtrList<surfaceScalarField> phiFluid;
  PtrList<uniformDimensionedVectorField> gFluid;
  PtrList<uniformDimensionedScalarField> hRefFluid;
  PtrList<volScalarField> ghFluid;
  PtrList<surfaceScalarField> ghfFluid;
  PtrList<compressible::turbulenceModel> turbulence;
  PtrList<volScalarField> p_rghFluid;
  PtrList<radiationModel> radiation;
  PtrList<volScalarField> KFluid;
  PtrList<volScalarField> dpdtFluid;
  List<scalar> initialMassFluid;
  List<bool> frozenFlowFluid;
  PtrList<IOMRFZoneList> MRFfluid;
  PtrList<fv::options> fluidFvOptions;
  // readPIMPLEControls
  int nOuterCorr;
  bool correctPhi;
  //. initContinuityErrs
  scalarList cumulativeContErr;
  // readTimeControls
  scalar DiNum;
  // compressibleMultiRegionCourantNo
  scalar CoNum;
  // createTimeControls
  bool adjustTimeStep;
  scalar maxCo;
  scalar maxDeltaT;

 public:

  //- Runtime type information
  TypeName("chtMultiRegionFoam");

  // Constructors

  //- Construct from mesh and phase name
  chtMultiRegionFoam(Time& runTime);

  //- Destructor
  virtual ~chtMultiRegionFoam();

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
#include "chtMultiRegionFoam.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
