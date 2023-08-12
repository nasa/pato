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

\*---------------------------------------------------------------------------*/

#include "elasticOrthoSolidFoamSolidMechanicsModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::elasticOrthoSolidFoamSolidMechanicsModel::elasticOrthoSolidFoamSolidMechanicsModel
(
    const fvMesh& mesh_foam,
    const word& regionName
)
  :
simpleSolidMechanicsModel(mesh_foam, regionName),
materialFailureCriteriaModel_autoPtr(simpleMaterialFailureCriteriaModel::New(mesh_foam, regionName)),
materialFailureCriteriaModel_ptr(&materialFailureCriteriaModel_autoPtr()),
materialFailureMassRemovalModel_autoPtr(simpleMaterialFailureMassRemovalModel::New(mesh_foam, regionName)),
materialFailureMassRemovalModel_ptr(&materialFailureMassRemovalModel_autoPtr()),
update_every_n_iter(materialDict_.subDict("SolidMechanics").lookupOrDefault<int>("update_every_n_iter",1)),
print_every_iCorr_iter(materialDict_.subDict("SolidMechanics").lookupOrDefault<int>("print_every_iCorr_iter",0)),
shear_stress_patch_names(materialDict_.subDict("SolidMechanics").lookupOrDefault<wordList>("shear_stress_patch_names",wordList())),
shear_stress_field_names(materialDict_.subDict("SolidMechanics").lookupOrDefault<wordList>("shear_stress_field_names",wordList())),
correct_pressure_BC(false),
iter(0),
dynamicMesh(isA<dynamicFvMesh>(mesh_foam)),
mesh_foam_(mesh_foam),
mesh_extend
(
    lookup_foam_extend_mesh(mesh_foam,regionName)
),
mesh(mesh_extend),
runTime(const_cast<Foam_extend_::Time&>(mesh.time())),
U
(
    Foam_extend_::IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::MUST_READ,
        Foam_extend_::IOobject::AUTO_WRITE
    ),
    mesh
),
gradU
(
    Foam_extend_::IOobject
    (
        "grad(U)",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::NO_READ,
        Foam_extend_::IOobject::NO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedTensor("zero", Foam_extend_::dimless, Foam_extend_::tensor::zero)
),
epsilon
(
    Foam_extend_::IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::READ_IF_PRESENT,
        Foam_extend_::IOobject::AUTO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedSymmTensor("zero", Foam_extend_::dimless, Foam_extend_::symmTensor::zero)
),
sigma(Foam_extend_::meshLookupOrConstructSymmTensor(mesh,"sigma",Foam_extend_::dimensionedSymmTensor("zero", Foam_extend_::dimForce/Foam_extend_::dimArea, Foam_extend_::symmTensor::zero))),
divSigmaExp
(
    Foam_extend_::IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::NO_READ,
        Foam_extend_::IOobject::NO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedVector("zero", Foam_extend_::dimForce/Foam_extend_::dimVolume, Foam_extend_::vector::zero)
),
rheology(sigma, U),
C(rheology.C()),
Cf(Foam_extend_::fvc::interpolate(C, "C")),
K(rheology.K()),
Kf(Foam_extend_::fvc::interpolate(K, "K")),
rho(rheology.rho()),
n(mesh.Sf()/mesh.magSf()),
divSigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaExp")),
massModel_(refModel<simpleMassModel>()),
p(massModel_.refVolField<scalar>("p"))
{
  {
#   include "using_extend.H"
#   include "elasticOrthoSolidFoam_createHistory.H"
#   include "elasticOrthoSolidFoam_readDivSigmaExpMethod.H"

    solidInterfaceCorr = rheology.solidInterfaceActive();

    solidInterfacePtr = nullptr;

    {
      if (solidInterfaceCorr) {
        // constitutiveModel is now in charge of solidInterface
        solidInterfacePtr = &rheology.solInterface();
        solidInterfacePtr->modifyProperties(Cf, Kf);

        //- solidInterface needs muf and lambdaf to be used for divSigmaExp
        if (divSigmaExpMethod != "surface") {
          FatalError
              << "divSigmaExp must be 'surface' when solidInterface is on"
              << exit(FatalError);
        }

        // check grad scheme
        if (word(mesh.schemesDict().gradSchemes().lookup("grad(U)"))
            != "leastSquaresSolidInterface") {
          Warning
              << "The grad(U) gradScheme should be "
              << "leastSquaresSolidInterface for the solidInterface "
              << "procedure to be correct"
              << endl;
        }
      }
    }
  }

  // Shear patch names
  forAll(shear_stress_patch_names, i) {
    int index = mesh_extend.boundaryMesh().findPatchID(shear_stress_patch_names[i]);
    if (index < 0) {
      FatalError << shear_stress_patch_names[i] << " not found in boundary mesh." << exit(FatalError);
    }
    if(!Foam_extend_::isA<Foam_extend_::solidTractionFvPatchVectorField>(U.boundaryField()[index])) {
      FatalError << "The type of \"" << shear_stress_patch_names[i] <<
                 "\" boundary for the field U has to be \"solidTractionFvPatchVectorField\"." << exit(FatalError);
    }
  }

  // Shear patch fields
  if (shear_stress_patch_names.size() != 0) {
    if (shear_stress_field_names.size() != 3) {
      FatalError << "shear_patch_names.size() != 3" << exit(FatalError);
    }
    shear_stress_fields.resize(shear_stress_field_names.size());
    forAll(shear_stress_field_names, i) {
      if (shear_stress_field_names[i] == "None") {
        word axis = "x";
        if (i==1) axis = "y";
        if (i==2) axis = "z";
        shear_stress_fields.set(i,
                                new volScalarField
                                (
                                    IOobject
                                    (
                                        "no_shear_stress_"+axis,
                                        mesh_foam.time().timeName(),
                                        mesh_foam,
                                        IOobject::NO_READ,
                                        IOobject::NO_WRITE
                                    ),
                                    mesh_foam,
                                    dimensionedScalar("0",dimMass/dimLength/pow(dimTime,2),0)
                                )
                               );
      } else {
        shear_stress_fields.set(i,
                                new volScalarField
                                (
                                    IOobject
                                    (
                                        shear_stress_field_names[i],
                                        mesh_foam.time().timeName(),
                                        mesh_foam,
                                        IOobject::MUST_READ,
                                        IOobject::NO_WRITE
                                    ),
                                    mesh_foam
                                )
                               );
      }
    }

    // Correct pressure BC
    correct_pressure_BC = true;
    wordList modelNames = {"Mass","Energy"};
    forAll(modelNames, i) {
      if (materialDict_.isDict(modelNames[i])) {
        const word modelType = materialDict_.subDict(modelNames[i]).lookupOrDefault<word>(modelNames[i]+"Type","no");
        if (modelType != "no") {
          correct_pressure_BC = false;
        }
      }
    }
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elasticOrthoSolidFoamSolidMechanicsModel::~elasticOrthoSolidFoamSolidMechanicsModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::elasticOrthoSolidFoamSolidMechanicsModel::update()
{
  if (iter%update_every_n_iter==0) {
    updateBC();
    iter=0;
#   include "using_extend.H"
#   include "elasticOrthoSolidFoam_readSolidMechanicsControls.H"

    int iCorr = 0;
    int iCorr_print = 0;
    lduSolverPerformance solverPerf;
    scalar initialResidual = 1.0;
    scalar relativeResidual = 1.0;
    lduMatrix::debug = 0;

    do {
      U.storePrevIter();

#     include "elasticOrthoSolidFoam_calculateDivSigmaExp.H"

      //- Linear momentum equation
      fvVectorMatrix UEqn
      (
          rho*Foam_extend_::fvm::d2dt2(U)
          ==
          Foam_extend_::fvm::laplacian(Kf, U, "laplacian(K,U)")
          + divSigmaExp
      );

      if (solidInterfaceCorr) {
        solidInterfacePtr->correct(UEqn);
      }

      solverPerf = UEqn.solve();

      if (iCorr == 0) {
        initialResidual = solverPerf.initialResidual();
      }

      U.relax();

      gradU = Foam_extend_::fvc::grad(U); // use leastSquaresSolidInterface

#     include "elasticOrthoSolidFoam_calculateRelativeResidual.H"

      if (iCorr % infoFrequency == 0) {
        Info<< "\tTime " << runTime.value()
            << ", Corr " << iCorr
            << ", Solving for " << U.name()
            << " using " << solverPerf.solverName()
            << ", res = " << solverPerf.initialResidual()
            << ", rel res = " << relativeResidual
            << ", inner iters " << solverPerf.nIterations() << endl;
      }

      if (print_every_iCorr_iter>0 && Pstream::nProcs() == 1) {
        if (iCorr_print%print_every_iCorr_iter==0) {
          Info << "\tPriting Corr " << iCorr << endl;
          iCorr_print=0;
          wordList fieldNames = {"U","sigma","interpolate(K)","divSigmaExp","rho","rheologyLawStoredC","K"};
          U.write();
          sigma.write();
          Kf.write();
          divSigmaExp.write();
          rho.write();
          C.write();
          K.write();
          fileName dirCorr = runTime.timeName()+"/"+regionName_+"/corr";
          if (!isDir(dirCorr)) {
            system("mkdir -p "+dirCorr);
          }
          forAll(fieldNames, i) {
            fileName old_path = runTime.timeName()+"/"+regionName_+"/"+fieldNames[i];
            fileName new_path = dirCorr + "/" + fieldNames[i] +"_corr_"+std::to_string(iCorr);
            if (isFile(old_path)) {
              word cmd = "mv " + old_path + " " + new_path;
              cmd.replaceAll("(","\\(").replaceAll(")","\\)");
              system(cmd.c_str());
            }
          }
        }
        iCorr_print++;
      }
    } while
    (
        solverPerf.initialResidual() > convergenceTolerance
        && ++iCorr < nCorr
    );

    Info<< nl << "Time " << runTime.value() << ", Solving for " << U.name()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations " << iCorr
        << endl;

#   include "elasticOrthoSolidFoam_calculateEpsilonSigma.H"
#   include "elasticOrthoSolidFoam_writeFields.H"
#   include "elasticOrthoSolidFoam_writeHistory.H"
    materialFailureCriteriaModel_ptr->updateFailureCriteria(); // update the failure criteria
    materialFailureMassRemovalModel_ptr->updateMassRemoval(); // update mass removal
  } else {
    if (iter%update_every_n_iter==1) {
      materialFailureCriteriaModel_ptr->cleanFields(); // clean the failure criteria field
      materialFailureMassRemovalModel_ptr->cleanFields(); // clean mass removal fields
    }
  }
  iter++;
}

void Foam::elasticOrthoSolidFoamSolidMechanicsModel::updateBC()
{
  // Update the shear stress BC
  forAll(shear_stress_fields, i) {
    shear_stress_fields[i].correctBoundaryConditions();
  }

  // Update the pressure BC
  if (correct_pressure_BC) {
    p.correctBoundaryConditions();
  }

  // Update U solidTraction BC
  forAll(shear_stress_patch_names, j) {
    int index = mesh_extend.boundaryMesh().findPatchID(shear_stress_patch_names[j]);
    Foam_extend_::solidTractionFvPatchVectorField& U_st=Foam_extend_::refCast<\
        Foam_extend_::solidTractionFvPatchVectorField>(U.boundaryField()[index]);

    // Update the pressure
    Foam_extend_::scalarField& pressure = U_st.pressure();
    forAll(pressure, faceI) {
      pressure[faceI] = p.boundaryField()[index][faceI];
    }

    // Update the traction
    Foam_extend_::vectorField& traction = U_st.traction();
    forAll(traction, faceI) {
      Foam_extend_::vector vec(0,0,0);
      vec.x() = shear_stress_fields[0].boundaryField()[index][faceI];
      vec.y() = shear_stress_fields[1].boundaryField()[index][faceI];
      vec.z() = shear_stress_fields[2].boundaryField()[index][faceI];
      traction[faceI] = vec;
    }
  }
}


// ************************************************************************* //
