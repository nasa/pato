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
    Foam::simpleEnergyModel

Description

SourceFiles
    simpleEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleEnergyModel_H
#define simpleEnergyModel_H

#include "simpleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class simpleEnergyModel
:
public simpleModel
{

 private:

  // Private Member Functions

 protected:

  // Protected data

 public:

  //- Runtime type information
  TypeName("simpleEnergyModel");

  //- Model name
  static const word modelName;

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleEnergyModel,
      fvMesh,
      (const fvMesh& mesh,const word& regionName),
      (mesh, regionName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleEnergyModel
  (
      const fvMesh& mesh,
      const word& regionName
  );


  //- Selector
  static autoPtr<simpleEnergyModel> New
  (
      const fvMesh& mesh,
      const word& regionName
  );

  //- Destructor
  virtual ~simpleEnergyModel();

  //- Loop over all cells and update all thermodynamic properties
  virtual void update() = 0;

};

#include "simpleEnergyModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
