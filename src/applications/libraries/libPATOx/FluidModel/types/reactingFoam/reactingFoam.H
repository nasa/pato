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
    Foam::reactingFoam

Description
    Tabulated thermo in function of pressure and temperature

SourceFiles
    reactingFoam.H
    reactingFoam.C

\*---------------------------------------------------------------------------*/

#ifndef reactingFoam_H
#define reactingFoam_H

#include "basicFluidModel.H"
// ADDED INCLUDE
#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IFstream.H"
#include "IOstream.H"
#include "ISstream.H"
#include "RectangularMatrix.H"
#include "fixedValueFvPatchField.H"
#include "fixedGradientFvPatchField.H"
//#include "psiReactionThermo.H"
#include "CombustionModel.H"
//#include "thermophysicalModels/chemistryModel/lnInclude/psiChemistryModel.H"
#include "thermophysicalModels/chemistryModel/lnInclude/chemistrySolver.H"
#include "OFstream.H"
#include "thermoPhysicsTypes.H"
#include "basicMultiComponentMixture.H"
#include "probes.H"
#include "IOprobes.H"
#undef Log // conflict between OpenFoam and Mutation++ on the alias "Log" -> undefined here from OpenFoam.
#include <mutation++/mutation++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "localEulerDdtScheme.H"
#include "coupledMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class reactingFoam
:
public basicFluidModel
{

 protected:
  // createFluidMeshes
  PtrList<dynamicFvMesh> fluidRegions;
  // createFluidFields
  PtrList<rhoReactionThermo> thermoFluid;
  PtrList<CombustionModel<rhoReactionThermo>> reactionFluid;
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
//  PtrList<volScalarField> dQFluid;
  List<scalar> initialMassFluid;
  List<bool> frozenFlowFluid;
  PtrList<IOMRFZoneList> MRFfluid;
  PtrList<fv::options> fluidFvOptions;
  PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable> fieldsFluid;
  // initContinuityErrs
  scalarList cumulativeContErr;
  // readTimeControls
  scalar DiNum;
  // compressibleMultiRegionCourantNo
  scalar CoNum;
  // createTimeControls
  bool adjustTimeStep;
  scalar maxCo;
  scalar maxDeltaT;
  // readPIMPLEControls
  int nOuterCorr;
  bool correctPhi;
  PtrList<volScalarField> diffY;

 public:

  //- Runtime type information
  TypeName("reactingFoam");

  // Constructors

  //- Construct from mesh and phase name
  reactingFoam(Time& runTime);

  //- Destructor
  virtual ~reactingFoam();

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
#include "reactingFoam.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
