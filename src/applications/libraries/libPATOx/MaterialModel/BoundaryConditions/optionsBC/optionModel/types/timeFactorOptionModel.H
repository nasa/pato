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
    Foam::timeFactorOptionModel

Description
    Foam::timeFactorOptionModel

    This class provides the factor option for boundary conditions.

SourceFile
    timeFactorOptionModel.C

\*---------------------------------------------------------------------------*/

#ifndef timeFactorOptionModel_H
#define timeFactorOptionModel_H

#include "simpleOptionModel.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class timeFactorOptionModel Declaration
\*---------------------------------------------------------------------------*/

class timeFactorOptionModel
:
public simpleOptionModel
{

  //- Time factor file
  fileName timeFactorFile_;

  //- List of the field name and the factor index
  List<Tuple2<word,label>> timeFactorList_;

  //- Time factor data
  List<scalarList> timeFactorData_;

 public:

  // Constructors

  //- Construct from patch and internal field
  timeFactorOptionModel(const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName);

  //- Runtime type information
  TypeName("timeFactorOption");

  // Member functions

  //- Update option
//  virtual void update();

  //- Update option for faceI
  virtual void update(label faceI);

  //- Print the fields and factors
  void print();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
