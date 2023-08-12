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
    Foam::DarcyForchheimerLawMassModel

Description
    Semi-implicit formulation of mass conservation (with Darcy-Forchheimer's law
-ie. , momentum conservation- substituted inside the mass conservation equation)

SourceFiles
    DarcyForchheimerLawMassModel.H
    DarcyForchheimerLawMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef DarcyForchheimerLawMassModel_H
#define DarcyForchheimerLawMassModel_H

#include "simpleMassModel.H"
#include "simpleEnergyModel.H"
#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DarcyForchheimerLawMassModel Declaration
\*---------------------------------------------------------------------------*/

class DarcyForchheimerLawMassModel
:
public simpleMassModel
{

 protected:

  //- Forchheimer's coefficient [1/m]
  volTensorField& Beta;

  //- Gravity field [m/s2]
  const uniformDimensionedVectorField& g;

  //- Solid permeability [m2]
  volTensorField& K;

  //- Gaseous mass flux [kg/m2/s]
  volVectorField& mDotG;

  //- Interpolation on the faces of the gas mass flux [kg/m3/s]
  surfaceVectorField& mDotGFace;

  //- Pyrolysis gas mass flow rate at the walls [kg/m2/s]
  volScalarField& mDotGw;

  //- Gaseous pressure [Pa]
  volScalarField& p;

  //- Gaseous velocity times gas volume fraction [m/s]
  volVectorField& U;

  //- Gaseous velocity [m/s]
  volVectorField& v;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Energy Model
  const simpleEnergyModel& energyModel_;

  //- Gaseous temperature [K]
  const volScalarField& T;

  //- GasPropertiesModel
  const simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gaseous volume fraction [-]
  const volScalarField& eps_g;

  //- Gaseous molar mass [kg/kmol]
  const volScalarField& M;

  //- Gaseous viscosity [kg/m/s]
  const volScalarField& mu;

  //- Gaseous mass density [-]
  volScalarField& rho_g;

  //- PyrolysisModel
  const simplePyrolysisModel& pyrolysisModel_;

  //- Total gas pyrolysis production [kg/m3/s]
  const volScalarField& piTotal;

  //- Part of the convective mass flux [s] : q_conv = dx . ( Gamma dx p ) + dx . (Gamma_G g)
  volTensorField& Gamma;

  //- Symmetrical of the part of the convective mass flux [s] : q_conv = dx . ( Gamma dx p ) + dx . (Gamma_G g)
  volSymmTensorField& Gamma_symm;

  //- Part of the storage mass flux [s] : q_sto = dt( Eta p )
  volScalarField& Eta;

  //- Part of the convective mass flux [(kg s)/m3] : q_conv = dx . ( Gamma dx p ) + dx . (Gamma_G g)
  volTensorField& Gamma_G;

  //- Interpolation on the face of the gas mass flux coming from gravity field [kg/m3/s]
  surfaceScalarField& phiG;

 public:

  //- Runtime type information
  TypeName("DarcyForchheimerLaw");

  // Constructors

  //- Construct from mesh and phase name
  DarcyForchheimerLawMassModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~DarcyForchheimerLawMassModel();

  //- Member Functions

  //- Update the mass conservation fields
  virtual void update();

  //- Update fields before to solve the mass conservation equation
  void beforeSolve();

  //- Update fields after to solve the mass conservation equation
  void afterSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
