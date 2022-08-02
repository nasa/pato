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
    Foam::rhoCentralFoamUserThermo

Description

SourceFiles
    rhoCentralFoamUserThermo.H
    rhoCentralFoamUserThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoCentralFoamUserThermo_H
#define rhoCentralFoamUserThermo_H

#include "basicFluidModel.H"
#include "fvCFD.H"
#include "basicUserThermo.H"
#include "turbulentUserThermo.H"
#include "fixedRhoFvPatchScalarField.H"
#include "rhoCentralFoam/directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlSolver Declaration
\*---------------------------------------------------------------------------*/

class rhoCentralFoamUserThermo
:
public basicFluidModel
{

 protected:

  //- fluid region
  wordList regionList;

  //- Mesh reference
  autoPtr<dynamicFvMesh> meshPtr;
  dynamicFvMesh& mesh;

  // Courant numbers used to adjust the time-step
  scalar CoNum;
  scalar meanCoNum;

  // createFields
//  autoPtr<psiThermo> pThermo;
  autoPtr<basicUserThermo> pThermo;
  basicUserThermo& thermo;
  volScalarField& e;
  volVectorField U;
  volScalarField rho;
  volVectorField rhoU;
  volScalarField rhoE;
  surfaceScalarField pos;
  surfaceScalarField neg;
  surfaceScalarField phi;
  autoPtr<compressibleUser::turbulenceModel> turbulence;

  //createFieldRefs
  volScalarField& p;
  const volScalarField& T;
  const volScalarField& psi;
  const volScalarField& mu;
  bool inviscid;

  // readFluxScheme
  word fluxScheme;

  // readTimeControls
  scalar DiNum;

  // createTimeControls
  bool adjustTimeStep;
  scalar maxCo;
  scalar maxDeltaT;

  // createRDeltaT
  bool LTS;
  tmp<volScalarField> trDeltaT;

  // rhoCentralFoam
  dimensionedScalar v_zero;
  surfaceScalarField rho_pos;
  surfaceScalarField rho_neg;
  surfaceVectorField rhoU_pos;
  surfaceVectorField rhoU_neg;
  surfaceVectorField U_pos;
  surfaceVectorField U_neg;
  surfaceScalarField aphiv_pos;
  surfaceScalarField aphiv_neg;
  surfaceScalarField phiv_pos;
  surfaceScalarField phiv_neg;
  volScalarField rPsi;
  surfaceScalarField rPsi_pos;
  surfaceScalarField rPsi_neg;
  volScalarField c;
  surfaceScalarField cSf_pos;
  surfaceScalarField cSf_neg;
  surfaceScalarField ap;
  surfaceScalarField am;
  surfaceScalarField a_pos;
  surfaceScalarField a_neg;
  surfaceScalarField p_pos;
  surfaceScalarField p_neg;
  surfaceScalarField e_pos;
  surfaceScalarField e_neg;
  surfaceScalarField aSf;


 public:

  //- Runtime type information
  TypeName("rhoCentralFoamUserThermo");

  // Constructors

  //- Construct from mesh and phase name
  rhoCentralFoamUserThermo(Time& runTime);

  //- Destructor
  virtual ~rhoCentralFoamUserThermo();

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
#include "rhoCentralFoamUserThermo.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
