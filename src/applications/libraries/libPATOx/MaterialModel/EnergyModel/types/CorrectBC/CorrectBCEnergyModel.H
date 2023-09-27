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
    along with OpenFOAM.  If CorrectBCt, see <http://www.gnu.org/licenses/>.

Class
    Foam::CorrectBCEnergyModel

Description

SourceFiles
    CorrectBCEnergyModel.H
    CorrectBCEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef CorrectBCEnergyModel_H
#define CorrectBCEnergyModel_H

#include "simpleEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CorrectBCEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class CorrectBCEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Temperature [K]
  volScalarField& T_;

 public:

  //- Runtime type information
  TypeName("CorrectBC");

  // Constructors

  //- Construct from mesh and phase name
  CorrectBCEnergyModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~CorrectBCEnergyModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
