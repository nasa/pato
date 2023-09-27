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
    Foam::simpleSolidMechanicsModel

Description

SourceFiles
    simpleSolidMechanicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleSolidMechanicsModel_H
#define simpleSolidMechanicsModel_H

#include "simpleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleSolidMechanicsModel Declaration
\*---------------------------------------------------------------------------*/

class simpleSolidMechanicsModel
:
public  simpleModel
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Plane stress ?
  Switch planeStress_;

 public:

  //- Runtime type information
  TypeName("simpleSolidMechanicsModel");

  //- Model name
  static const word modelName;

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleSolidMechanicsModel,
      fvMesh,
      (const fvMesh& mesh,const word& regionName),
      (mesh, regionName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleSolidMechanicsModel
  (
      const fvMesh& mesh,
      const word& regionName
  );


  //- Selector
  static autoPtr<simpleSolidMechanicsModel> New
  (
      const fvMesh& mesh,
      const word& regionName
  );

  //- Destructor
  virtual ~simpleSolidMechanicsModel();

  //- Update
  virtual void update() = 0;

  //- Thermal expansion coefficient [1/K]
  inline volScalarField& alpha() const;

  //- Jomaa Function (pyrolysis shrinkage) [-]
  inline volScalarField& xi() const;

  //- Plane stress ?
  inline Switch planeStress() const;


};

#include "simpleSolidMechanicsModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
