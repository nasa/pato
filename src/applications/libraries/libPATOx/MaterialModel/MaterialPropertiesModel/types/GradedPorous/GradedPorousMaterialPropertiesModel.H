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
    along with OpenFOAM.  If Porous_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::GradedPorousMaterialPropertiesModel

Description

SourceFiles
    GradedPorousMaterialPropertiesModel.H
    GradedPorousMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef GradedPorousMaterialPropertiesModel_H
#define GradedPorousMaterialPropertiesModel_H

#include "simpleMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "LinearInterpolationMaterialPropertiesObject.H"
#include "dimensionedTensor.H"
#include "IOFunctions.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GradedPorousMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class GradedPorousMaterialPropertiesModel
:
public simpleMaterialPropertiesModel
{

 protected:

  //- Material Properties directory
  const fileName materialPropertiesDirectory;

  //- Dictionary constantProperties
  IOdictionary constantPropertiesDictionary;

  //- Number of solid phases
  scalar& nSolidPhases_;

  //- Initial solid phase volume fractions [-]
  PtrList<volScalarField>& solidEpsI_;

  //- Initial solid phase densities [kg/m3]
  PtrList<volScalarField>& solidRhoI_;

  //- Solid phase densities [kg/m3]
  PtrList<volScalarField>& solidRho_;

  //- Solid phase volume fractions [-]
  PtrList<volScalarField>& solidEps_;

  //- File names of the grading properties (distant, volume fraction)
  PtrList<fileName> gradingFiles_;

  //- Patch names for the grading
  PtrList<word> gradingPatchNames_;

  //- Patch ID for the grading
  PtrList<label> gradingPatchID_;

  //- Switch between virgin or char material properties
  const word virginOrChar;

  //- True if there are several regions
  bool multiRegions;

  //- Total solid density [kg/m3]
  volScalarField& rho_s_;

  //- Virgin object (linear interpolation)
  autoPtr<LinearInterpolationMaterialPropertiesObject> virginObject;

  //- Char object (linear interpolation)
  autoPtr<LinearInterpolationMaterialPropertiesObject> charObject;

  //- Permeability [m2]
  volTensorField& K_;

  //- Absorptivity [-]
  volScalarField& absorptivity_;

  //- Emissivity [-]
  volScalarField& emissivity_;

  //- Temperature [K]
  volScalarField& T_;

  //- Pressure [Pa]
  volScalarField& p_;

  //- Char enthalpy (used in the surface mass and energy balance) [J/kg]
  volScalarField& h_c_;

  //- Averaged solid enthalpy (used in the energy model) [J/kg]
  volScalarField& h_bar_;

  //- Solid heat capacity [J/kg/K]
  volScalarField& cp_;

  //- Solid thermal conductivity [W/m/K]
  volTensorField& k_;

  //- Thermal conductivity in basis ijk (main directions of the conductivity tensor) [W/m/K]
  volTensorField kijk_;

  //- Symmetric thermal conductivity for multi-material [W/m/K]
  volSymmTensorField k_abl_sym_;

  //- Thermal conductivity coefficient in I direction
  dimensionedScalar kiCoef_;

  //- Thermal conductivity coefficient in J direction
  dimensionedScalar kjCoef_;

  //- Thermal conductivity coefficient in K direction
  dimensionedScalar kkCoef_;

  //- Tensor for the passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  dimensionedTensor tP_;

  //- Char permeability [m2]
  dimensionedTensor Kc_;

  //- Virgin permeability [m2]
  dimensionedTensor Kv_;

  //- Pyrolysis energy [W/m3]
  volScalarField& pyrolysisFlux_;

  //- detailedSolidEnthalpies flag (default false)
  Switch detailedSolidEnthalpies_;

  //- Reference to Pyrolysis Model
  simplePyrolysisModel& PyrolysisModel_;

  //- Pyrolysis reaction rates per reactions [kg/m3/s]
  PtrList<volScalarField>& piPyroReac_;

  //- Totalpyrolysis reaction rates [kg/m3/s]
  volScalarField& piTotal_;

  //- Advancement of pyrolysis reactions
  volScalarField& tau_;

  //- Pyrolysis reaction enthalpies [J/kg]
  PtrList<dimensionedScalar>& hp_;

  //- Virgin solid density [kg/m3]
  volScalarField& rho_v_;

  //- Char solid density [kg/m3]
  volScalarField& rho_c_;

  //- Initialization switch
  const Switch init_;

 public:

  //- Runtime type information
  TypeName("GradedPorous");

  // Constructors

  //- Construct from mesh and phase name
  GradedPorousMaterialPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~GradedPorousMaterialPropertiesModel();

  //- Member Functions

  //- Initialize the model
  Switch init();

  //- Read initial solid phase volume fractions
  void readEpsI();

  //- Read initial solid phase densities
  void readRhoI();

  //- Read char or virgin objects (linear interpolation)
  LinearInterpolationMaterialPropertiesObject readCharredVirginDict(word charredVirgin);

  //- Update emissivity and absorptivity
  void updateEmissivityAbsorptivity();

  //- Update the porous material properties
  virtual void update();

  //- Update char enthalpy (used in Bprime boundary conditions)
  void updateHc();

  //- Update averaged solid enthalpy (used in energy model)
  void updateHbar();

  //- Update solid heat capacity
  void updateCp();

  //- Update solid thermal conductivity
  void updatek();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
