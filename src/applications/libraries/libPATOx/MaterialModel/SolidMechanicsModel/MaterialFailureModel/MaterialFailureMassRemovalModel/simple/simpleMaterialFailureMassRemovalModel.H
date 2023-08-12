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
    Foam::simpleMaterialFailureMassRemovalModel

Description

SourceFiles
    simpleMaterialFailureMassRemovalModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMaterialFailureMassRemovalModel_H
#define simpleMaterialFailureMassRemovalModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleMaterialFailureMassRemovalModel Declaration
\*---------------------------------------------------------------------------*/

class simpleMaterialFailureMassRemovalModel
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

 public:

  //- Runtime type information
  TypeName("simpleMaterialFailureMassRemovalModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMaterialFailureMassRemovalModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMaterialFailureMassRemovalModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simpleMaterialFailureMassRemovalModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simpleMaterialFailureMassRemovalModel();

  //- Update material failure mass removal
  virtual void updateMassRemoval() = 0;

  //- Clean the fields
  virtual void cleanFields() = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
