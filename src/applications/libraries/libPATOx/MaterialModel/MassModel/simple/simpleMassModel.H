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
    Foam::simpleMassModel

Description

SourceFiles
    simpleMassModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMassModel_H
#define simpleMassModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleMassModel Declaration
\*---------------------------------------------------------------------------*/

class simpleMassModel
:
public  IOdictionary
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

  //- Dynamic mesh flag
  const bool dynamicMesh_;

  //- Pyrolysis gaseous mass flow rate [kg/m2/s]
  volVectorField& mDotG_;

  //- Gaseous velocity [m/s]
  volVectorField& vG_;

  //- Interpolation on the faces of the pyrolysis gaseous mass flow rate [kg/m2/s]
  surfaceVectorField& mDotGFace_;

 public:

  //- Runtime type information
  TypeName("simpleMassModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMassModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMassModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simpleMassModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simpleMassModel();

  //- Update
  virtual void update() = 0;

  //- Gaseous mass flux [kg/m2/s]
  inline volVectorField& mDotG() const;

  //- Gaseous velocity [m/s]
  inline volVectorField& vG() const;

  //- Interpolation on the faces of the gaseous mass flux
  inline surfaceVectorField& mDotGFace() const;

};

#include "simpleMassModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
