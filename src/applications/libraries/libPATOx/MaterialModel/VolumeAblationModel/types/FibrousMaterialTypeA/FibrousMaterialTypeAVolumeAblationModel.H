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
    Foam::FibrousMaterialTypeAVolumeAblationModel

Description

SourceFiles
    FibrousMaterialTypeAVolumeAblationModel.H
    FibrousMaterialTypeAVolumeAblationModel.C

\*---------------------------------------------------------------------------*/

#ifndef FibrousMaterialTypeAVolumeAblationModel_H
#define FibrousMaterialTypeAVolumeAblationModel_H

#include "simpleVolumeAblationModel.H"
#include "IOFunctions.H"
#include "simplePyrolysisModel.H"
#include "simpleMaterialChemistryModel.H"
#include "simpleGasPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FibrousMaterialTypeAVolumeAblationModel Declaration
\*---------------------------------------------------------------------------*/

class FibrousMaterialTypeAVolumeAblationModel
:
public simpleVolumeAblationModel
{
 protected:

  //- energy conservation type
  const word energyConservation_;

  //- Specific surface [1/m]
  volScalarField& specificSurface;

  //- Site density [mol/m2]
  volScalarField& siteDensity;

  //- Initial site density [mol/m2]
  volScalarField& siteDensityIni;

  //- Fiber radius [m]
  volScalarField& rT;

  //- Initial fiber radius [m]
  volScalarField& rT0;

  //- Solid volume fraction [-]
  volScalarField& eps_s_;

  //- Initial solid volume fraction [-]
  volScalarField& eps_s0_;

  //- Charred ablation rate [kg/s]
  dimensionedScalar m_dot_ablation;

  //- Initial fiber radius [m]
  const dimensionedScalar& rf0;

  //- Final fiber radius [m]
  const dimensionedScalar& rff;

  //- Failure fiber radius [m]
  const dimensionedScalar& rfFail;

  //- Fiber phase
  label fiberPhase;

  //- Matrix phase
  label matrixPhase;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Solid phase volume fractions [-]
  PList<volScalarField> solidEps_;

  //- Initial solid phase volume fractions [-]
  PList<volScalarField> solidEpsI_;

  //- Solid phase densities [kg/m3]
  PList<volScalarField> solidRho_;

  //- Initial solid phase densities [kg/m3]
  PList<volScalarField> solidRhoI_;

  //- External solid density [kg/m3]
  volScalarField rho_ext;

  //- MaterialChemistry Model
  const simpleMaterialChemistryModel& MaterialChemistryModel_;

  //- Heterogeneous rate [kg/m3/s]
  const volScalarField& omegaHeterogeneousRate_;

  //- Mass Model
  const simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gas density [kg/m3]
  const volScalarField& rho_g_;

  //- Molar mass of solid carbon [kg/kmol]
  const dimensionedScalar M_Cs_;

 public:

  //- Runtime type information
  TypeName("FibrousMaterialTypeA");

  //- Construct from mesh and phase name
  FibrousMaterialTypeAVolumeAblationModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~FibrousMaterialTypeAVolumeAblationModel();

  // Member Functions

  //- Update the volume ablation fields
  virtual void update();

  //- Update the solid carbon mass fraction
  virtual void updateSolidCarbonMassFraction(volScalarField& YiCs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
