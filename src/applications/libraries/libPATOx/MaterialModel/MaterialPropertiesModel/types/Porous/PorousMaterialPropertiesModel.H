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
    Foam::PorousMaterialPropertiesModel

Description

SourceFiles
    PorousMaterialPropertiesModel.H
    PorousMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef PorousMaterialPropertiesModel_H
#define PorousMaterialPropertiesModel_H

#include "simpleMaterialPropertiesModel.H"
#include "CommonMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "LinearInterpolationMaterialPropertiesObject.H"
#include "dimensionedTensor.H"
#include "IOFunctions.H"
#include "PMap.H"
#include "MapField.H"

#define bindFun(fun,name) std::bind(&Foam::PorousMaterialPropertiesModel::fun, this, name)
#define newMapField(name,type,init,update) new MapField(name,type,bindFun(init,name),bindFun(update,name))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PorousMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class PorousMaterialPropertiesModel
:
public simpleMaterialPropertiesModel,
public CommonMaterialPropertiesModel
{

 protected:

  // Protected Member Variables

  //- Flag for initVirginCharFiles function
  Switch initVirginCharFiles_;

  //- Virgin object (linear interpolation)
  autoPtr<LinearInterpolationMaterialPropertiesObject> virginObject;

  //- Char object (linear interpolation)
  autoPtr<LinearInterpolationMaterialPropertiesObject> charObject;

  //- Mass Model
  const simpleMassModel& massModel_;

  //- Pressure [Pa]
  const volScalarField& p_;

  //- Energy Model
  const simpleEnergyModel& energyModel_;

  //- Temperature [K]
  const volScalarField& T_;

  //- Solid mass density [kg/m3]
  volScalarField& rho_s_;

  //- Pyrolysis Model
  const simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of pyrolysis reactions [-]
  const volScalarField& tau_;

  //- Solid virgin mass density [kg/m3]
  const volScalarField& rho_v_;

  //- Solid charred mass density [kg/m3]
  const volScalarField& rho_c_;

  //- True if there are several regions
  bool multiRegions;

  //- Solid phase densities [kg/m3]
  PList<volScalarField> solidRho_;

  //- Solid phase volume fractions [-]
  PList<volScalarField> solidEps_;

  //- Pyrolysis gas production per reaction [kg/m3/s]
  PList<volScalarField> piPyroReac_;

 public:

  //- Runtime type information
  TypeName("Porous");

  // Constructors

  //- Construct from mesh and phase name
  PorousMaterialPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PorousMaterialPropertiesModel();

  // Member Functions

  // Initialization functions

  //- Init the porous material properties
  virtual void initialize();

  //- Init solid density [kg/m3]
  void init_rho_s(word name);

  //- Init fields using virgin/char files
  void initVirginCharFiles(word name);

  //- Init scalar properties
  void initScalarConstant(word name);

  //- Init scalar properties for linear interpolation (virgin/char)
  void initScalarLinear(word name);

  //- Init vector properties for linear interpolation (virgin/char)
  void initVectorLinear(word name);

  //- Init tensor properties for linear interpolation (virgin/char)
  void initTensorLinear(word word);

  //- Read char or virgin objects (linear interpolation)
  LinearInterpolationMaterialPropertiesObject readCharredVirginDict(word charredVirgin);

  // Update functions

  //- Update the porous material properties
  virtual void update();

  //- Update scalar properties
  void updateScalarConstant(word name);

  //- Init scalar properties for linear interpolation (virgin/char)
  void updateScalarLinear(word name);

  //- Init vector properties for linear interpolation (virgin/char)
  void updateVectorLinear(word name);

  //- Init tensor properties for linear interpolation (virgin/char)
  void updateTensorLinear(word name);

  //- Update emissivity and absorptivity
  void updateEmissivityAbsorptivity(word name);

  //- Update char enthalpy (used in Bprime boundary conditions)
  void updateHc(word name);

  //- Update averaged solid enthalpy and pyrolysis flux (used in energy model)
  void updatePyrolysisFlux(word name);

  //- Update solid density
  void update_rho_s(word name);

  //- Update solid heat capacity
  void update_cp(word name);

  //- Update solid thermal conductivity
  void updatek(word name);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
