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
    along with OpenFOAM.  If Porous_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::MapField

Description

SourceFiles
    MapField.H
    MapField.C

\*---------------------------------------------------------------------------*/

#ifndef MapField_H
#define MapField_H

#include "fvCFD.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MapField Declaration
\*---------------------------------------------------------------------------*/

class MapField
{

 public:

  // Constructors

  //- Construct from fieldName, fieldTypeName, f_init and f_update
  MapField(word fieldName,word fieldTypeName,
           std::function<void()> f_init,
           std::function<void()> f_update);

  // Public Members

  //- Field Name
  word fieldName_;

  //- Field Type Name
  word fieldTypeName_;

  //- Init function
  std::function<void()> f_init_;

  //- Update function
  std::function<void()> f_update_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
