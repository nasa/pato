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
    Foam::elasticSolidFoamSolidMechanicsModel

Description

SourceFiles
    elasticSolidFoamSolidMechanicsModel.H
    elasticSolidFoamSolidMechanicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef elasticSolidFoamSolidMechanicsModel_H
#define elasticSolidFoamSolidMechanicsModel_H

#define namespaceFoam
#include "simpleMaterialsModel.H"
#include "simpleSolidMechanicsModel.H"
#undef FALSE
#undef TRUE
#define namespaceFoam_extend_
#include "declare_macros_foam_extend.H" // Redeclare Foam Extend 4.1 macros
#include "Foam_extend_fvCFD.H"
#include "Foam_extend_dynamicFvMesh.H"
#include "Foam_extend_constitutiveModel.H"
#include "Foam_extend_solidInterface.H"
#include "declare_macros_openfoam.H" // Redeclare OpenFOAM 7 macros

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class elasticSolidFoamSolidMechanicsModel Declaration
\*---------------------------------------------------------------------------*/

class elasticSolidFoamSolidMechanicsModel
:
public simpleSolidMechanicsModel
{

 protected:

  //- Update every N iterations
  int update_every_n_iter;

  //- Iteration number
  int iter;

  // Foam_extend_ variables
  Switch dynamicMesh;
  Foam_extend_::fvMesh& mesh_extend;
  Foam_extend_::fvMesh& mesh;
  Foam_extend_::Time& runTime;
  Foam_extend_::uniformDimensionedVectorField g;
  Foam_extend_::volVectorField U;
  Foam_extend_::volTensorField gradU;
  Foam_extend_::surfaceVectorField snGradU;
  Foam_extend_::volVectorField V;
  Foam_extend_::volTensorField gradV;
  Foam_extend_::surfaceVectorField snGradV;
  Foam_extend_::volSymmTensorField epsilon;
  Foam_extend_::volSymmTensorField sigma;
  Foam_extend_::volVectorField divSigmaExp;
  Foam_extend_::constitutiveModel rheology;
  Foam_extend_::volScalarField rho;
  Foam_extend_::volScalarField mu;
  Foam_extend_::volScalarField lambda;
  Foam_extend_::surfaceScalarField muf;
  Foam_extend_::surfaceScalarField lambdaf;
  Foam_extend_::surfaceVectorField n;
  Foam_extend_::volVectorField aitkenDelta;
  scalar aitkenInitialRes;
  scalar aitkenTheta;
  Foam_extend_::OFstream* filePtr;
  label historyPatchID;
  word historyPatchName;
  Foam_extend_::word divSigmaExpMethod;
  bool solidInterfaceCorr;
  Foam_extend_::solidInterface* solidInterfacePtr;

 public:

  //- Runtime type information
  TypeName("elasticSolidFoam");

  // Constructors

  //- Construct from mesh and phase name
  elasticSolidFoamSolidMechanicsModel(const fvMesh& mesh_foam,const word& regionName);

  //- Destructor
  virtual ~elasticSolidFoamSolidMechanicsModel();

  //- Member Functions

  //- Update
  virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
