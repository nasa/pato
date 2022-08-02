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
    Foam::DarcyLaw2TMassModel

Description
    Semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
    momentum conservation- substituted inside the mass conservation equation)

SourceFiles
    DarcyLaw2TMassModel.H
    DarcyLaw2TMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef DarcyLaw2TMassModel_H
#define DarcyLaw2TMassModel_H

#include "simpleMassModel.H"
#include "simpleGasPropertiesModel.H"
#include "simpleMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DarcyLaw2TMassModel Declaration
\*---------------------------------------------------------------------------*/

class DarcyLaw2TMassModel
:
public simpleMassModel
{

 protected:

  //- Mesh reference
  const fvMesh& mesh_;

  //- Dictionary name
  const word& dictName_;

  //- GasPropertiesModel
  simpleGasPropertiesModel& gasPropertiesModel_;

  //- MaterialPropertiesModel
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- PyrolysisModel
  simplePyrolysisModel& pyrolysisModel_;

  //- Gaseous pressure [Pa]
  volScalarField& p;

  //- Gaseous temperature [K]
  volScalarField& Tg_;

  //- Gaseous molar mass [kg/kmol]
  volScalarField& M;

  //- Gaseous viscosity [kg/m/s]
  volScalarField& mu;

  //- Gaseous volume fraction [-]
  volScalarField& eps_g;

  //- Gaseous mass density [-]
  volScalarField& rho_g;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Solid permeability [m2]
  volTensorField& K;

  //- Total gas pyrolysis production [kg/m3/s]
  volScalarField& piTotal;

  //- Part of the convective mass flux [s] : q_conv = dx . ( Gamma dx p )
  volTensorField Gamma;

  //- Symmetrical of the part of the convective mass flux [s] : q_conv = dx . ( Gamma dx p )
  volSymmTensorField Gamma_symm;

  //- Part of the storage mass flux [s] : q_sto = dt( Beta p )
  volScalarField Beta;

  //- Gaseous velocity [m/s]
  volVectorField& v;

  //- Gaseous velocity times gas volume fraction [m/s]
  volVectorField& U;

  //- Gaseous mass flux [kg/m2/s]
  volVectorField& mDotG;

  //- Interpolation on the faces of the gas mass flux [kg/m3/s]
  surfaceVectorField& mDotGFace;

 public:

  //- Runtime type information
  TypeName("DarcyLaw2T");

  // Constructors

  //- Construct from mesh and phase name
  DarcyLaw2TMassModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~DarcyLaw2TMassModel();

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
