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
    along with OpenFOAM.  If Pyrolysist, see <http://www.gnu.org/licenses/>.

Class
    Foam::Pyrolysis2TEnergyModel

Description

SourceFiles
    Pyrolysis2TEnergyModel.H
    Pyrolysis2TEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef Pyrolysis2TEnergyModel_H
#define Pyrolysis2TEnergyModel_H

#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "simpleMaterialPropertiesModel.H"
#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Pyrolysis2TEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class Pyrolysis2TEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Reference to the region name
  const word& dictName_;

  //- Mass Model
  simpleMassModel& MassModel_;

  //- Material Properties Model
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- Material Properties directory
  const fileName materialPropertiesDirectory;

  //- Dictionary constantProperties
  IOdictionary constantPropertiesDictionary;

  //- Gas Properties Model
  simpleGasPropertiesModel& gasPropertiesModel_;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Gas enthalpy [J/kg]
  volScalarField& h_g;

  //- Pressure [Pa]
  volScalarField& p;

  //- Molar Mass [kg/kmol]
  volScalarField& M;

  //- Viscosity [kg/m/s]
  volScalarField& mu;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Gas temperature [K]
  volScalarField& Tg;

  //- Solid temperature [K]
  volScalarField& Ts;

  //- Permeability [m2]
  const volTensorField& K;

  //- Gas volume fraction (1 - porosity) [-]
  const volScalarField& eps_g;

  //- Gas density [kg/m3]
  const volScalarField& rho_g;

  //- Solid density [kg/m3]
  volScalarField& rho_s;

  //- Gas heat capacity [J/kg/K]
  volScalarField& cp_g;

  //- Solid heat capacity [J/kg/K]
  volScalarField& cp_s;

  //- Pyrolysis flux [W/m3] - updated in Material Properties Model
  volScalarField& pyrolysisFlux_;

  //- Gaseous velocity [m/s]
  volVectorField& vg;

  //- Gaseous thermal conductivity [W/m/K]
  volScalarField& k_g;

  //- Equivalent thermal conductivity [W/m/K]
  volTensorField& k_a;

  //- Solid thermal conductivity [W/m/K]
  volTensorField k_s;

  //- Gas volume fraction [-] * Gas density [kg/m3] * Gas enthalpy [J/kg]
  volScalarField epsgRhog;

  //- Convection term / grad(pressure)
  volTensorField GammaHg;

  //- Gaseous mass flux [kg/m2/s]
  surfaceScalarField phi_g;

  //- Convective heat exchange between solid and gas phases
  dimensionedScalar Hv0_;
  //const volScalarField& Hv0_;

 public:

  //- Runtime type information
  TypeName("Pyrolysis2T");

  // Constructors

  //- Construct from mesh and phase name
  Pyrolysis2TEnergyModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~Pyrolysis2TEnergyModel();

  //- Member Functions

  //- Update the energy equation fields
  virtual void update();

  //- Before to solve the energy equation
  void beforeSolve();

  //- After to solve the energy equation
  void afterSolve();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
