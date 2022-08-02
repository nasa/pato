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
    Foam::simpleMaterialPropertiesModel

Description

SourceFiles
    simpleMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMaterialPropertiesModel_H
#define simpleMaterialPropertiesModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class simpleMaterialPropertiesModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Dictionary name
  const word& dictName_;

  //- Material dictionary
  IOdictionary materialDict_;

  //- Debug Switch
  const Switch debug_;

  //- dynamic mesh flag
  const bool dynamicMesh_;

  //- solid heat capacity [J/kg/K]
  volScalarField& cp_;

  //- solid thermal conductivity [W/m/K]
  volTensorField& k_;

  //- Solid densities [kg/m^3]
  PtrList<volScalarField> solidRho_;

  //- Solid volume fraction [-]
  PtrList<volScalarField> solidEps_;

  //- Solid initial volume fraction [-]
  PtrList<volScalarField> solidEpsI_;

  //- Solid initial densities [kg/m^3]
  PtrList<volScalarField> solidRhoI_;

  //- Permeability [m2]
  volTensorField& K_;

  //- Solid density [kg/m^3]
  volScalarField& rho_s_;

  //- Averaged solid enthalpy [J/kg]
  volScalarField& h_bar_;

  //- Emissivity [-]
  volScalarField& emissivity_;

  //- Absorptivity [-]
  volScalarField& absorptivity_;

  //- charred enthalpy [J/kg]
  volScalarField& h_c_;

  //- Pyrolysis flux [W/m3]
  volScalarField& pyrolysisFlux_;

  //- number of solid phases
  scalar nSolidPhases_;

 public:

  //- Runtime type information
  TypeName("simpleMaterialPropertiesModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMaterialPropertiesModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMaterialPropertiesModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simpleMaterialPropertiesModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simpleMaterialPropertiesModel();

  //- Update
  virtual void update() = 0;

  //- solid heat capacity [J/kg/K]
  inline volScalarField& cp() const;

  //- solid thermal conductivity [W/m/K]
  inline volTensorField& k() const;

  //- Solid densities [kg/m^3]
  inline PtrList<volScalarField>& solidRho();

  //- Solid volume fraction [-]
  inline PtrList<volScalarField>& solidEps();

  //- Solid initial volume fraction [-]
  inline PtrList<volScalarField>& solidEpsI() const;

  //- Solid initial densities [kg/m^3]
  inline PtrList<volScalarField>& solidRhoI() const;

  //- Permeability [m2]
  inline volTensorField& K() const;

  //- Solid density [kg/m^3]
  inline volScalarField& rho_s() const;

  //- Averaged solid enthalpy [J/kg]
  inline volScalarField& h_bar() const;

  //- Emissivity [-]
  inline volScalarField& emissivity() const;

  //- Absorptivity [-]
  inline volScalarField& absorptivity() const;

  //- charred enthalpy [J/kg]
  inline volScalarField& h_c() const;

  //- Pyrolysis flux [W/m3]
  inline volScalarField& pyrolysisFlux() const;

  //- number of solid phases
  inline scalar& nSolidPhases() const;

};

#include "simpleMaterialPropertiesModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
