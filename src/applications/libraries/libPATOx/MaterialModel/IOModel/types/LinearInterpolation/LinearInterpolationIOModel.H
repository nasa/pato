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
    along with OpenFOAM.  If LinearInterpolationt, see <http://www.gnu.org/licenses/>.

Class
    Foam::LinearInterpolationIOModel

Description

SourceFiles
    LinearInterpolationIOModel.H
    LinearInterpolationIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef LinearInterpolationIOModel_H
#define LinearInterpolationIOModel_H

#include "noIOModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class LinearInterpolationIOModel Declaration
\*---------------------------------------------------------------------------*/

class LinearInterpolationIOModel
:
public noIOModel
{

 protected:

  //- Start model initialization after parent class
  const Switch startModelInit_;

  //- Output word list
  const wordList outputList_;

  //- setLinear switch
  const Switch setLinear_;

  //- List of the linear field names
  const wordList linearFieldList_;

  //- Initilialization flag
  const bool initOutput_;

 public:

  //- Runtime type information
  TypeName("LinearInterpolation");

  // Constructors

  //- Construct from mesh and phase name
  LinearInterpolationIOModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~LinearInterpolationIOModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
