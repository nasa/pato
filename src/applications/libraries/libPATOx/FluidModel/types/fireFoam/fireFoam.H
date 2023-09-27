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
    Foam::fireFoam

Description

SourceFiles
    fireFoam.H
    fireFoam.C

\*---------------------------------------------------------------------------*/

#ifndef fireFoam_H
#define fireFoam_H

#include "basicFluidModel.H"
// ADDED INCLUDE
#include "turbulentFluidThermoModel.H"
#include "basicReactingCloud.H"
#include "surfaceFilmModel.H"
//#include "pyrolysisModelCollection.H"
#include "radiationModel.H"
#include "SLGThermo.H"
//#include "solidChemistryModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fvOptions.H"
#include "compressibleCourantNo.H"
#include "coordinateSystem.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "coupledMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class fireFoam
:
public basicFluidModel
{

 protected:
  // createFluidMeshes
  PtrList<dynamicFvMesh> fluidRegions;
  // createFluidFields
  PtrList<psiReactionThermo> pThermo;
  PtrList<volScalarField> rhoFluid;
  PtrList<volVectorField> UFluid;
  PtrList<surfaceScalarField> phiFluid;
  PtrList<uniformDimensionedVectorField> gFluid;
  PtrList<uniformDimensionedScalarField> hRefFluid;
  PtrList<uniformDimensionedScalarField> pRefFluid;
  PtrList<volScalarField> ghFluid;
  PtrList<surfaceScalarField> ghfFluid;
  PtrList<compressible::turbulenceModel> turbulence;
  PtrList<CombustionModel<psiReactionThermo>> combustion;
  PtrList<volScalarField> p_rghFluid;
  PtrList<volScalarField> ph_rghFluid;
  PtrList<radiationModel> radiation;
  PtrList<volScalarField> KFluid;
  PtrList<volScalarField> dpdtFluid;
  PtrList<IOMRFZoneList> MRFfluid;
  PtrList<fv::options> fluidFvOptions;
  PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable> fieldsFluid;
  PtrList<volScalarField> diffY;
  List<scalar> initialMassFluid;
  // readPIMPLEControls
  int nOuterCorr;
  int nHydrostaticCorrectors;
  bool correctPhi;
  bool hydrostaticInitialization;
  // initContinuityErrs
  scalarList cumulativeContErr;
  // readTimeControls
  scalar maxDi;
  scalar DiNum;
  // compressibleMultiRegionCourantNo
  scalar CoNum;
  // createTimeControls
  bool adjustTimeStep;
  scalar adjustStartTime;
  scalar maxDeltaT;
  scalar minDeltaT;
  scalar maxCo;

 public:

  //- Runtime type information
  TypeName("fireFoam");

  // Constructors

  //- Construct from mesh and phase name
  fireFoam(Time& runTime);

  //- Destructord
  virtual ~fireFoam();

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
#include "fireFoam.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
