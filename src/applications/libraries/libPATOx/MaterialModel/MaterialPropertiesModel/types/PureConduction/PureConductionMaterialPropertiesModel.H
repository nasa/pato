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
    along with OpenFOAM.  If PureConduction_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::PureConductionMaterialPropertiesModel

Description
    PureConduction material properties (k, cp ,rho_s) with inteprolation from virgin and char

SourceFiles
    PureConductionMaterialPropertiesModel.H
    PureConductionMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef PureConductionMaterialPropertiesModel_H
#define PureConductionMaterialPropertiesModel_H

#include "simpleMaterialPropertiesModel.H"
#include "LinearInterpolationMaterialPropertiesObject.H"
#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PureConductionMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class PureConductionMaterialPropertiesModel
:
public simpleMaterialPropertiesModel
{

 protected:

  //- Material Properties directory
  const fileName materialPropertiesDirectory;

  //- Dictionary constantProperties
  IOdictionary constantPropertiesDictionary;

  //- Switch between virgin or char material properties
  const word virginOrChar;

  //- Number of solid phases
  scalar& nSolidPhases_;

  //- Initial solid phase volume fractions [-]
  PtrList<volScalarField>& solidEpsI_;

  //- Initial solid phase densities [kg/m3]
  PtrList<volScalarField>& solidRhoI_;

  //- Virgin object (linear interpolation)
  autoPtr<LinearInterpolationMaterialPropertiesObject> virginObject;

  //- Char object (linear interpolation)
  autoPtr<LinearInterpolationMaterialPropertiesObject> charObject;

  //- Thermal conductivity coefficient in I direction
  dimensionedScalar kiCoef_;

  //- Thermal conductivity coefficient in J direction
  dimensionedScalar kjCoef_;

  //- Thermal conductivity coefficient in K direction
  dimensionedScalar kkCoef_;

  //- Tensor for the passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  dimensionedTensor tP_;

  //- Solid phase densities [kg/m3]
  PtrList<volScalarField>& solidRho_;

  //- Solid phase volume fractions [-]
  PtrList<volScalarField>& solidEps_;

  //- Reference to Pyrolysis Model
  simplePyrolysisModel& PyrolysisModel_;

  //- Advancement of pyrolysis reactions
  volScalarField& tau_;

  //- Total solid density [kg/m3]
  volScalarField& rho_s_;

  //- Virgin solid density [kg/m3]
  volScalarField& rho_v_;

  //- Char solid density [kg/m3]
  volScalarField& rho_c_;

  //- Temperature [K]
  volScalarField& T_;

  //- Pressure [Pa]
  volScalarField& p_;

  //- Solid heat capacity [J/kg/K]
  volScalarField& cp_;

  //- Solid thermal conductivity [W/m/K]
  volTensorField& k_;

  //- Absorptivity [-]
  volScalarField& absorptivity_;

  //- Emissivity [-]
  volScalarField& emissivity_;

  //- Char enthalpy (used in the surface mass and energy balance) [J/kg]
  volScalarField& h_c_;

  //- Thermal conductivity in basis ijk (main directions of the conductivity tensor) [W/m/K]
  volTensorField kijk_;

  //- Symmetric thermal conductivity for multi-material [W/m/K]
  volSymmTensorField k_abl_sym_;

 public:

  //- Runtime type information
  TypeName("PureConduction");

  // Constructors

  //- Construct from mesh and phase name
  PureConductionMaterialPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PureConductionMaterialPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Read initial solid volume fraction
  inline void readEpsI();

  //- Read initial solid mass density
  inline void readRhoI();

  //- Read the charred and virgin dictionaries
  LinearInterpolationMaterialPropertiesObject readCharredVirginDict(word charredVirgin);

  //- Update the heat capacity
  inline void updateCp();

  //- Update the thermal conductivity
  inline void updatek();

  //- Update the emissivity and absorptivity
  inline void updateEmissivityAbsorptivity();

  //- Update the charred enthalpy
  inline void updateHc();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
