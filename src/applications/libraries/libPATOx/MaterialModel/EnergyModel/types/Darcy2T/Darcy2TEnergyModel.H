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
    Foam::Darcy2TEnergyModel

Description

SourceFiles
    Darcy2TEnergyModel.H
    Darcy2TEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef Darcy2TEnergyModel_H
#define Darcy2TEnergyModel_H

#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "simpleGasPropertiesModel.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Darcy2TEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class Darcy2TEnergyModel
:
public simpleEnergyModel
{

 protected:

  // - Source term (OpenFOAM fvOptions)
  fv::options& fvOptions_;

  //- Solid heat capacity [J/kg/K]
  const volScalarField& cp_s;

  //- Convective heat exchange between solid and gas phases
  volScalarField&  Hv_;

  //- Equivalent thermal conductivity [W/m/K]
  volTensorField& k_a;

  //- Pyrolysis flux [W/m3] - updated in Material Properties Model
  const volScalarField& pyrolysisFlux_;

  //- Solid density [kg/m3]
  const volScalarField& rho_s;

  //- Gas temperature [K]
  volScalarField& Tg;

  //- Solid temperature [K]
  volScalarField& Ts;

  //- Particle size for sphere's bed
  const dimensionedScalar& Dp_;

  //- Constant convective heat exchange between solid and gas phases
  const dimensionedScalar& Hv0_;

  //- Heat transfer model type
  const word HvType_;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Mass Model
  const simpleMassModel& massModel_;

  //- Gravity field [m/s2]
  const uniformDimensionedVectorField g;

  //- Permeability [m2]
  const volTensorField& K;

  //- Pressure [Pa]
  const volScalarField& p;

  //- Gaseous velocity [m/s]
  const volVectorField& vg;

  //- Gas Properties Model
  const simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gas heat capacity [J/kg/K]
  const volScalarField& cp_g;

  //- Gas volume fraction (1 - porosity) [-]
  const volScalarField& eps_g;

  //- Gas enthalpy [J/kg]
  const volScalarField& h_g;

  //- Gaseous thermal conductivity [W/m/K]
  const volScalarField& k_g;

  //- Molar Mass [kg/kmol]
  const volScalarField& M;

  //- Viscosity [kg/m/s]
  const volScalarField& mu;

  //- Gas density [kg/m3]
  const volScalarField& rho_g;

  //- Gas volume fraction [-] * Gas density [kg/m3] * Gas enthalpy [J/kg]
  volScalarField& epsgRhog;

  //- Convection term (gravitational part)
  volTensorField& Gamma_GHg;

  //- Solid thermal conductivity [W/m/K]
  volTensorField& k_s;

  //- Interpolation on the face of the convection term coming from gravity field [kg/m3/s]
  surfaceScalarField& phiHg;

  //- Gaseous mass flux [kg/m2/s]
  surfaceScalarField& phi_g;

 public:

  //- Runtime type information
  TypeName("Darcy2T");

  // Constructors

  //- Construct from mesh and phase name
  Darcy2TEnergyModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~Darcy2TEnergyModel();

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
