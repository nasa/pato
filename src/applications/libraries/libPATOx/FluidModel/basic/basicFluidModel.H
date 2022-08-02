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
    Foam::basicFluidModel

Description

SourceFiles
    basicFluidModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicFluidModel_H
#define basicFluidModel_H

#include "fvCFD.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class basicFluidModel Declaration
\*---------------------------------------------------------------------------*/

class basicFluidModel
{

 private:

  // Private Member Functions



 protected:
  //- Region properties
  regionProperties rp;

  //- runTime reference
  Time& runTime;

  //- Delay solid computation
  const scalar delaySolid_;

  //- Delay solid computation
  const scalar delayFluid_;

  //- Criteria of stopping solid update and beginning fluid update
  const scalar solidStitchTolerance_;

  //- Criteria of stopping fluid update and beginning solid update
  const scalar fluidStitchTolerance_;

  //- Absolute Criteria of stopping solid update and beginning fluid update
  const scalar forceFluidUpdate_;

  //- Fluid stitch criteria
  scalar maxDeltaTw_;

 public:

  //- Runtime type information
  TypeName("basicFluidModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      basicFluidModel,
      Time,
      (Time& runTime),
      (runTime)
  );


  // Constructors

  //- Construct from mesh and phase name
  basicFluidModel
  (
      Time& runTime
  );


  //- Selector
  static autoPtr<basicFluidModel> New
  (
      Time& runTime
  );

  //- Destructor
  virtual ~basicFluidModel();

  //- Update before runTime++ in the while loop
  virtual void updateBefore() = 0;

  //- Update after runTime++ in the while loop
  virtual void updateAfter() = 0;

  //- Delay solid computation
  inline const scalar& delaySolid();

  //- Delay Fluid computation
  inline const scalar& delayFluid();

  //- Criteria of stopping solid update and beginning fluid update
  inline const scalar& solidStitchTolerance();

  //- Criteria of stopping fluid update and beginning solid update
  inline const scalar& fluidStitchTolerance();

  //- Absolute Criteria of stopping solid update and beginning fluid update
  inline const scalar& forceFluidUpdate();

  //- Fluid Stitch criteria
  inline scalar& maxDeltaTw();

  //- Used to initialize/read the thermodynamic variables

  //- Creates a scalar field with MUST_READ if not found in the mesh database
  inline volScalarField& meshLookupOrConstructScalar
  (
      const fvMesh& mesh,
      const char* name
  );

  //- Creates a scalar field with READ_IF_PRESENT if not found in the mesh database
  inline volScalarField& meshLookupOrConstructScalar
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedScalar& dimScalar_
  );

  //- Creates a scalar field with READ_IF_PRESENT and a given boundary condition for all patches if not found in the mesh database
  inline volScalarField& meshLookupOrConstructScalar
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedScalar& dimScalar_,
      const wordList& BC_
  );

  //- Creates a vector field with MUST_READ if not found in the mesh database
  inline volVectorField& meshLookupOrConstructVector
  (
      const fvMesh& mesh,
      const char* name
  );

  //- Creates a vector field  with READ_IF_PRESENT for all patches if not found in the mesh database
  inline volVectorField& meshLookupOrConstructVector
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedVector& dimVec_
  );

  //- Creates a vector field with READ_IF_PRESENT and a given boundary condition for all patches if not found in the mesh database
  inline volVectorField& meshLookupOrConstructVector
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedVector& dimVector_,
      const wordList& BC_
  );

  //- Creates a tensor field if not found in the mesh database
  inline volTensorField& meshLookupOrConstructTensor
  (
      const fvMesh& mesh,
      const char* name
  );

  //- Creates a tensor field with given a boundary condition for all patches if not found in the mesh database
  inline volTensorField& meshLookupOrConstructTensor
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedTensor& dimTens_
  );

  //- Creates a tensor field with READ_IF_PRESENT and a given boundary condition for all patches if not found in the mesh database
  inline volTensorField& meshLookupOrConstructTensor
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedTensor& dimTensor_,
      const wordList& BC_
  );

};

#include "basicFluidModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
