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
    Foam::simplePyrolysisModel

Description

SourceFiles
    simplePyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef simplePyrolysisModel_H
#define simplePyrolysisModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simplePyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class simplePyrolysisModel
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

  //- Return the pyrolysis reaction rates per reactions [kg/m^3/s]
  PtrList<volScalarField> piPyroReac_;

  //- Pyrolysis reaction rates per specie/element [kg/m^3/s]
  PtrList<volScalarField> pi_;

  //- Total pyrolysis reaction rates [kg/m^3/s]
  volScalarField& piTotal_;

  //- Advancement of the pyrolysis reaction [-]
  volScalarField& tau_;

  //- Temperature [K]
  volScalarField& Ta_;

  //- Solid density [kg/m^3]
  volScalarField& rho_s_;

  //- Virgin solid density [kg/m^3]
  volScalarField& rho_v_;

  // Charred solid density [kg/m^3]
  volScalarField& rho_c_;

  //- Pyrolysis reaction enthalpies [J/kg] - used when detailedSolidEnthalpy == "yes" in Material Properties Model
  PtrList<dimensionedScalar> hp_;

 public:

  //- Runtime type information
  TypeName("simplePyrolysisModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simplePyrolysisModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simplePyrolysisModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simplePyrolysisModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simplePyrolysisModel();

  //- Update
  virtual void update() = 0;

  //- Return const access to the mesh database
  inline const fvMesh& mesh() const;

  //- Return the pyrolysis reaction rates per reactions [kg/m^3/s]
  inline PtrList<volScalarField>& piPyroReac();

  //- Return the pyrolysis reaction rates per specie/element [kg/m^3/s]
  inline PtrList<volScalarField>& pi();

  //- Return the total pyrolysis reaction rates [kg/m^3/s]
  inline volScalarField& piTotal();

  //- Advancement of the pyrolysis reaction [-]
  inline volScalarField& tau();

  //- Temperature [K]
  inline volScalarField& Ta();

  //- Solid density [kg/m^3]
  inline volScalarField& rho_s();

  //- Virgin solid density [kg/m^3]
  inline volScalarField& rho_v();

  //- Charred solid density [kg/m^3]
  inline volScalarField& rho_c();

  //- Pyrolysis reaction enthalpies [J/kg]
  inline PtrList<dimensionedScalar>& hp();

};

#include "simplePyrolysisModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
