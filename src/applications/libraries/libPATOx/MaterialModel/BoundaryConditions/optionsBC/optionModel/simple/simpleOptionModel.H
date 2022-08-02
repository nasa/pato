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
    Foam::simpleOptionModel

Description
    Foam::simpleOptionModel

    This class provides boundary condition options.

SourceFile
    simpleOptionModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleOptionModel_H
#define simpleOptionModel_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class simpleOptionModel Declaration
\*---------------------------------------------------------------------------*/

class simpleOptionModel
:
public IOdictionary
{
 protected:
  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Dictionary
  const dictionary dict_;

  //- Debug switch
  const label patchID_;

  //- Type name
  const word optionTypeName_;

 public:

  // Constructors

  //- Runtime type information
  TypeName("simpleOption");

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleOptionModel,
      fvMesh,
      (const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName),
      (mesh, dict, patchID, optionTypeName)
  );

  //- Construct from patch and internal field
  simpleOptionModel(const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName);

  //- Selector
  static autoPtr<simpleOptionModel> New
  (
      const fvMesh& mesh,
      const dictionary& dict,
      const label& patchID,
      const word& optionTypeName
  );

  // Member functions

  //- Update option
  virtual void update() {
    Info << "update simple" << endl;
  };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
