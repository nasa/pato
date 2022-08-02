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
    Foam::simpleBoundaryMappingModel

Description

SourceFiles
    simpleBoundaryMappingModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleBoundaryMappingModel_H
#define simpleBoundaryMappingModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"
#include "mathFunctions.H"
#include "fluxFactor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleBoundaryMappingModel Declaration
\*---------------------------------------------------------------------------*/

class simpleBoundaryMappingModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Fields needed
  const wordList neededFields_;

  //- Fields found in the dictionary as constant values
  List<Tuple2<word,scalar> > constantFields_;

  //- Dictionary name
  const dictionary dict_;

  //- Switch debug
  Switch debug_;

  //- Mapping file name
  const fileName mappingFileName_;

  //- Mapping information ( (fieldName index) ... )
  const List<Tuple2<word,word> > mappingFields_;

  //- Mapping fields name
  wordList mappingFieldsName_;

  //- Mapping fields column in mappingFileName_ data
  scalarList mappingFieldsColumn_;

  //- Current time (s) for all the mesh patches and data fields
  List<scalarList> currentTimePatchesDataFields_;

 public:

  //- Flux factor object
  fluxFactor * fluxFactor_ptr;

  //- Runtime type information
  TypeName("simpleBoundaryMappingModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleBoundaryMappingModel,
      fvMesh,
      (const fvMesh& mesh, const wordList& neededFields, const dictionary dict),
      (mesh, neededFields, dict)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleBoundaryMappingModel
  (
      const fvMesh& mesh,
      const wordList& neededFields,
      const dictionary dict
  );


  //- Selector
  static autoPtr<simpleBoundaryMappingModel> New
  (
      const fvMesh& mesh,
      const wordList& neededFields,
      const dictionary dict
  );

  //- Destructor
  virtual ~simpleBoundaryMappingModel();

  //- Update
  virtual void update(scalar timeValue, label patchID, word fieldName) = 0;

  //- write in os
  virtual void write(Ostream& os) const = 0;

  //- Initialize the boundary mapping model
  inline void init();

  //- Return mapping fields name
  inline wordList& mappingFieldsName() const;

};

#include "simpleBoundaryMappingModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
