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
    Foam::factorOptionModel

Description
    Foam::factorOptionModel

    This class provides the factor option for boundary conditions.

SourceFile
    factorOptionModel.C

\*---------------------------------------------------------------------------*/

#ifndef factorOptionModel_H
#define factorOptionModel_H

#include "simpleOptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class factorOptionModel Declaration
\*---------------------------------------------------------------------------*/

class factorOptionModel
:
public simpleOptionModel
{

  //- Current time
  scalar currentTime_;

//- List of the field name and the factor value
  List<Tuple2<word,scalar>> factorList;

 public:

  // Constructors

  //- Construct from patch and internal field
  factorOptionModel(const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName);

  //- Runtime type information
  TypeName("factorOption");

  // Member functions

  //- Update option
  virtual void update();

  //- Print the fields and factors
  void print();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
