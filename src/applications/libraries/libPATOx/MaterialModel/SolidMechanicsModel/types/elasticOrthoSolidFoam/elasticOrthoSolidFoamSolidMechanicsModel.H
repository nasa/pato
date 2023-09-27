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
    Foam::elasticOrthoSolidFoamSolidMechanicsModel

Description

SourceFiles
    elasticOrthoSolidFoamSolidMechanicsModel.H
    elasticOrthoSolidFoamSolidMechanicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef elasticOrthoSolidFoamSolidMechanicsModel_H
#define elasticOrthoSolidFoamSolidMechanicsModel_H

#define namespaceFoam
#include "simpleMaterialsModel.H"
#include "simpleSolidMechanicsModel.H"
#include "simpleMassModel.H"
#include "boundaryMappingFvPatchScalarField.H"
#include "IOFunctions.H"
#include "simpleMaterialFailureCriteriaModel.H"
#include "simpleMaterialFailureMassRemovalModel.H"
#include "IOFunctions_extend.H"
#undef FALSE
#undef TRUE
#define namespaceFoam_extend_
#include "declare_macros_foam_extend.H" // Redeclare Foam Extend 4.1 macros
#include "Foam_extend_fvCFD.H"
#include "Foam_extend_dynamicFvMesh.H"
#include "Foam_extend_constitutiveModel.H"
#include "Foam_extend_solidInterface.H"
#include "Foam_extend_solidTractionFvPatchVectorField.H"
#include "declare_macros_openfoam.H" // Redeclare OpenFOAM 7 macros

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class elasticOrthoSolidFoamSolidMechanicsModel Declaration
\*---------------------------------------------------------------------------*/

class elasticOrthoSolidFoamSolidMechanicsModel
:
public simpleSolidMechanicsModel
{

 protected:

  //- Material failure criteria model autoPtr
  autoPtr<simpleMaterialFailureCriteriaModel> materialFailureCriteriaModel_autoPtr;

  //- Material failure criteria model pointer
  simpleMaterialFailureCriteriaModel * materialFailureCriteriaModel_ptr;

  //- Material failure mass removal model autoPtr
  autoPtr<simpleMaterialFailureMassRemovalModel> materialFailureMassRemovalModel_autoPtr;

  //- Material failure mass removal model pointer
  simpleMaterialFailureMassRemovalModel * materialFailureMassRemovalModel_ptr;

  //- Update every N iterations
  const int update_every_n_iter;

  //- Print foam-extend fields every iCorr iterations
  const int print_every_iCorr_iter;

  //- List of the shear patch names
  const wordList shear_stress_patch_names;

  //- List of the shear stress field names
  const wordList shear_stress_field_names;

  //- List of the shear stress fields
  PtrList<volScalarField> shear_stress_fields;

  //- Flag for correct pressure bc
  Switch correct_pressure_BC;

  //- Iteration number
  int iter;

  //- Flag for fynamic mesh
  Switch dynamicMesh;

  //- Mesh from OpenFOAM
  const fvMesh& mesh_foam_;

  // Foam_extend_ variables
  Foam_extend_::fvMesh& mesh_extend;
  Foam_extend_::fvMesh& mesh;
  Foam_extend_::Time& runTime;
  Foam_extend_::volVectorField U;
  Foam_extend_::volTensorField gradU;
  Foam_extend_::volSymmTensorField epsilon;
  Foam_extend_::volSymmTensorField& sigma;
  Foam_extend_::volVectorField divSigmaExp;
  Foam_extend_::constitutiveModel rheology;
  Foam_extend_::volSymmTensor4thOrderField C;
  Foam_extend_::surfaceSymmTensor4thOrderField Cf;
  Foam_extend_::volDiagTensorField K;
  Foam_extend_::surfaceDiagTensorField Kf;
  Foam_extend_::volScalarField rho;
  Foam_extend_::surfaceVectorField n;
  Foam_extend_::word divSigmaExpMethod;
  bool solidInterfaceCorr;
  Foam_extend_::solidInterface* solidInterfacePtr;
  Foam_extend_::OFstream* filePtr;
  Foam_extend_::word historyPatchName;
  Foam_extend_::label historyPatchID;

  // Mass Model
  simpleMassModel& massModel_;

  // Pressure [Pa]
  volScalarField& p;

 public:

  //- Runtime type information
  TypeName("elasticOrthoSolidFoam");

  // Constructors

  //- Construct from mesh and phase name
  elasticOrthoSolidFoamSolidMechanicsModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~elasticOrthoSolidFoamSolidMechanicsModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Update the BC
  void updateBC();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
