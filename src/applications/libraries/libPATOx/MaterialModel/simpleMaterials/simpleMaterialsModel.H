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

#define namespaceFoam
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "regionProperties.H"
#include "IOFunctions.H"
#include "simpleMaterialModel.H"
#if defined(FOAM_EXTEND)
#undef FALSE
#undef TRUE
#define namespaceFoam_extend_
#include "declare_macros_foam_extend.H" // Redeclare Foam Extend 4.1 macros
#include "Foam_extend_fvCFD.H"
#include "Foam_extend_dynamicFvMesh.H"
#include "declare_macros_openfoam.H" // Redeclare OpenFOAM 7 macros
#endif
using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class simpleMaterialModel;

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

  //- Dynamic mesh list
  PtrList<dynamicFvMesh> dynamicMeshesList_;

#if defined(FOAM_EXTEND)
  //- Foam extend time
  Foam_extend_::Time foam_extend_runTime_;

  //- Foam extend meshes list
  PtrList<Foam_extend_::fvMesh> foam_extend_meshes_list_;
#endif

  // solid region names list
  wordList solidRegionNames_;

  //- Solid stitch criteria NAmrofel 26/11/2019
  scalar simpleMinDeltaTw_;

  //- Write fields flag in controlDict
  const Switch writeFields_;

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
  scalar& minDeltaTw();

  //- Write fields flag in controlDict
  const Switch& writeFields() const;

  //- Return meshes list
  const PtrList<fvMesh>& meshesList() const;

  //- Return dynamic meshes list
  const PtrList<dynamicFvMesh>& dynamicMeshesList() const;

#if defined(FOAM_EXTEND)
  //- Return foam extend time
  const Foam_extend_::Time& foam_extend_runTime() const;

  //- Return foam extend meshes list
  const PtrList<Foam_extend_::fvMesh>& foam_extend_meshes_list() const;

  //- Get the foam extend mesh for a given region
  const Foam_extend_::fvMesh& get_foam_extend_mesh(word regionName) const;
#endif

};

#include "simpleMaterialsModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
