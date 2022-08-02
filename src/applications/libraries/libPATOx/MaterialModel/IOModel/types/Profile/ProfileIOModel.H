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
    along with OpenFOAM.  If Profilet, see <http://www.gnu.org/licenses/>.

Class
    Foam::ProfileIOModel

Description

SourceFiles
    ProfileIOModel.H
    ProfileIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef ProfileIOModel_H
#define ProfileIOModel_H

#include "noIOModel.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ProfileIOModel Declaration
\*---------------------------------------------------------------------------*/

class ProfileIOModel
:
public noIOModel
{

 protected:

  const List<word> plot1DProfileList_;
  const Switch plot1DMassLoss_;
  const word plot1DFolder_;
  volScalarField& rho_s_;
  scalar V0_;
  scalar rho_s0_;
  const bool initOutput_;
  const Switch createFolder_;
  OFstream os_massLoss_;
  word topPatch_;
  word bottomPatch_;

 public:

  //- Runtime type information
  TypeName("Profile");

  // Constructors

  //- Construct from mesh and phase name
  ProfileIOModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~ProfileIOModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Create a folder
  Switch createFolder();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
