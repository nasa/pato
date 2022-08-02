/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::simpleGasPropertiesModel

Description

SourceFiles
    simpleGasPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleGasPropertiesModel_H
#define simpleGasPropertiesModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleGasPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class simpleGasPropertiesModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Material dictionary
  IOdictionary materialDict_;

  //- Debug Switch
  const Switch debug_;

  //- dynamic mesh flag
  const bool dynamicMesh_;

  //- Average molar mass of the gas phase [kg/kmol]
  volScalarField& M_;

  //- Average of gas heat capacity [J/kg/K]
  volScalarField& cp_g_;

  //- Average viscosity of the gas phase [kg/m/s]
  volScalarField& mu_;

  //- Average enthalpy of the gas phase [J/kg]
  volScalarField& h_g_;

  //- gas density [kg/m3]
  volScalarField& rho_g_;

  //- gas thermal conductivity [W/m-K]
  volScalarField& k_g_;

  //- gas volume fraction [-]
  volScalarField& eps_g_;

  //- Diffusion coefficients [m2/s]
  PtrList<volScalarField> Dm_;

  //- Flux transported by diffusion of the species [W/m3]
  volScalarField& Ediff_;

 public:

  //- Runtime type information
  TypeName("simpleGasPropertiesModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleGasPropertiesModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleGasPropertiesModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simpleGasPropertiesModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simpleGasPropertiesModel();

  //- Update
  virtual void update() = 0;

  //- Average molar mass of the gas phase [kg/kmol]
  inline volScalarField& M();

  //- Average of gas heat capacity [J/kg/K]
  inline volScalarField& cp_g();

  //- Average viscosity of the gas phase [kg/m/s]
  inline volScalarField& mu();

  //- Average enthalpy of the gas phase [J/kg]
  inline volScalarField& h_g();

  //- Density of the gas phase [kg/m3]
  inline volScalarField& rho_g();

  //- Thermal conductivity of the gas phase [W/m-K]
  inline volScalarField& k_g();

  //- Volume fraction of the gas phase [-]
  inline volScalarField& eps_g();

  //- Diffusion coefficients [m2/s]
  inline PtrList<volScalarField>& Dm();

  //- Flux transported by diffusion of the species [W/m3]
  inline volScalarField& Ediff();


};

#include "simpleGasPropertiesModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
