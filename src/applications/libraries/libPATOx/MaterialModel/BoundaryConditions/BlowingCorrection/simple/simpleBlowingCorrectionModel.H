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
    Foam::simpleBlowingCorrectionModel

Description

SourceFiles
    simpleBlowingCorrectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleBlowingCorrectionModel_H
#define simpleBlowingCorrectionModel_H

#include "fvCFD.H"
#include "simpleModel.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleBlowingCorrectionModel Declaration
\*---------------------------------------------------------------------------*/

class simpleBlowingCorrectionModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary name
  const word regionName_;

  //- Current patch ID
  const label currentPatchID_;

  //- Dictionary name
  const dictionary dict_;

  //- Initialized flag
  Switch initialized_;

  //- volScalarFields
  PMap<volScalarField> scalarFields_;

  //- volVectorFields
  PMap<volVectorField> vectorFields_;

  //- volTensorFields
  PMap<volTensorField> tensorFields_;

  //- Blowing correction type name
  const word blowingCorrectionType_;

  //- Debug flag
  const Switch debug_;

 public:

  //- Runtime type information
  TypeName("simpleBlowingCorrectionModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleBlowingCorrectionModel,
      fvMesh,
      (const fvMesh& mesh, const word& phaseName, const word& regionName, const label& currentPatchID, const dictionary dict),
      (mesh, phaseName, regionName, currentPatchID, dict)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleBlowingCorrectionModel
  (
      const fvMesh& mesh,
      const word& phaseName,
      const word& regionName,
      const label& currentPatchID,
      const dictionary dict
  );


  //- Selector
  static autoPtr<simpleBlowingCorrectionModel> New
  (
      const fvMesh& mesh,
      const word& phaseName,
      const word& regionName,
      const label& currentPatchID,
      const dictionary dict
  );

  //- Destructor
  virtual ~simpleBlowingCorrectionModel();

  //- Initialize
  virtual void initialize() = 0;

  //- Update
  virtual void update() = 0;

  //- write in os
  virtual void write(Ostream& os) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
