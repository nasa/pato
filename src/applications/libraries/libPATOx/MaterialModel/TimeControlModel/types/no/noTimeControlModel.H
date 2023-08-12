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
    Foam::noTimeControlModel

Description

SourceFiles
    noTimeControlModel.H
    noTimeControlModel.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeControlModel_H
#define noTimeControlModel_H

#include "simpleTimeControlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noTimeControlModel Declaration
\*---------------------------------------------------------------------------*/

class noTimeControlModel
:
public simpleTimeControlModel
{

 protected:
 public:

  //- Runtime type information
  TypeName("no");

  // Constructors

  //- Construct from mesh and phase name
  noTimeControlModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~noTimeControlModel();

  //- Member Functions

  //- Update thermodynamic and transport properties in the field
  virtual void update();

  //- Update stitch criteria
  virtual scalar updateMinDeltaTw();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
