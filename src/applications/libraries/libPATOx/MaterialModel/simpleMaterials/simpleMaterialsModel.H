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
    Foam::simpleMaterialsModel

Description

SourceFiles
    simpleMaterialsModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMaterialsModel_H
#define simpleMaterialsModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "regionProperties.H"
#include "IOFunctions.H"
#include "simpleMaterialModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleMaterialsModel Declaration
\*---------------------------------------------------------------------------*/

class simpleMaterialsModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Material Model list
  PtrList<simpleMaterialModel> materialsList_;

  //- Reference to the time
  const Time& runTime_;

  //- Region properties
  regionProperties rp;

  //- Mesh list
  PtrList<fvMesh> meshesList_;

  //- Mesh list
  PtrList<dynamicFvMesh> dynamicMeshesList_;

  // solid region names list
  wordList solidRegionNames_;

  //- Solid stitch criteria NAmrofel 26/11/2019
  scalar simpleMinDeltaTw_;

 public:

  //- Runtime type information
  TypeName("simpleMaterialsModel");

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMaterialsModel,
      Time,
      (const Time& runTime),
      (runTime)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMaterialsModel
  (
      const Time& runTime
  );

  //- Destructor
  virtual ~simpleMaterialsModel();

  //- Update all the materials
  void update();

  //- Fluid Stitch criteria, NAmrofel 26/11/2019
  inline scalar& minDeltaTw();

};

#include "simpleMaterialsModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
