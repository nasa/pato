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
    Foam::optionsBC

Description
    Foam::optionsBC

    This boundary condition provides common boundary condition options.

SourceFile
    optionsBC.C

\*---------------------------------------------------------------------------*/

#ifndef optionsBC_H
#define optionsBC_H

#include "fvCFD.H"
#include "simpleOptionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class optionsBC Declaration
\*---------------------------------------------------------------------------*/

class optionsBC
{
  //- List of the options (ptr)
  List<simpleOptionModel*> options_list_ptr;

  //- List of the options (autoPtr)
  List<autoPtr<simpleOptionModel> > options_list_autoPtr;

 public:

  // Constructors

  //- Construct from patch and internal field
  optionsBC(const fvMesh& mesh, const dictionary& dict, const label& patchID);

  // Member functions

  //- Update options
  void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
