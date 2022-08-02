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
    along with OpenFOAM.  If SpeciesConservationt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SpeciesConservationMaterialChemistryModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SpeciesConservationMaterialChemistryModel::SpeciesConservationMaterialChemistryModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialChemistryModel(mesh, dictName),
mixtureMutation(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").lookup("mixture")),
mix_(simpleMaterialChemistryModel::mixture_),
omegaHeterogeneousRate_(meshLookupOrConstructScalar(mesh, "omegaHeterogeneousRate", dimensionedScalar("0", dimMass/pow3(dimLength)/dimTime, 0.0))),
omegaHeterogeneousEnergy_(meshLookupOrConstructScalar(mesh, "omegaHeterogeneousEnergy", dimensionedScalar("0", dimMass/dimLength/pow3(dimTime), 0.0))),
materialPropertiesDirectory(fileName(simpleMaterialChemistryModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh.time().db().parent(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
elementNames_(simpleMaterialChemistryModel::elementNames_),
speciesNames_(simpleMaterialChemistryModel::speciesNames_),
massFractions_(simpleMaterialChemistryModel::massFractions_),
oldMassFractions_(simpleMaterialChemistryModel::oldMassFractions_),
moleFractions_(simpleMaterialChemistryModel::moleFractions_),
Ydefault_
(
    IOobject
    (
        "Ydefault",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh_
),
eta0_(constantPropertiesDictionary.lookup("eta0")),
T_(meshLookupOrConstructScalar(mesh,"Ta")),
p_(meshLookupOrConstructScalar(mesh,"p")),
h_ew_(meshLookupOrConstructScalar(mesh, "h_ew",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0)))),
pThermo_(psiReactionThermo::New(mesh)),
thermo_(pThermo_()),
pChemistry_(BasicFiniteRateChemistryModel<psiReactionThermo>::New(thermo_)),
chemistry_(pChemistry_()),
dtChem_(min(refCast<const BasicFiniteRateChemistryModel<psiReactionThermo>>(chemistry_).deltaTChem()[0],mesh.time().deltaT().value())),
//  thermo_.validate(args.executable(), "h");
composition_(thermo_.composition()),
Y_(composition_.Y()),
speciesIndexMutation_(simpleMaterialChemistryModel::speciesIndexMutation_),
init_(init()),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
pi_(pyrolysisModel_.pi()),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,dictName,"MaterialProperties")),
h_c_(materialPropertiesModel_.h_c()),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
rho_g_(gasPropertiesModel_.rho_g()),
eps_g_(gasPropertiesModel_.eps_g()),
Dm_(gasPropertiesModel_.Dm()),
massModel_(meshLookupOrConstructModel<simpleMassModel>(mesh,dictName,"Mass")),
mDotG_(massModel_.mDotG()),
mDotGFace_(massModel_.mDotGFace()),
mDotGw_(meshLookupOrConstructScalar(mesh, "mDotGw",dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, scalar(0.0)))),
volumeAblationModel_(meshLookupOrConstructModel<simpleVolumeAblationModel>(mesh,dictName,"VolumeAblation")),
volumeAblationType_(simpleMaterialChemistryModel::materialDict_.isDict("VolumeAblation")?word(simpleMaterialChemistryModel::materialDict_.subDict("VolumeAblation").lookup("VolumeAblationType")):""),
gasPropertiesType_(simpleMaterialChemistryModel::materialDict_.isDict("GasProperties")?word(simpleMaterialChemistryModel::materialDict_.subDict("GasProperties").lookup("GasPropertiesType")):"")
{

  word materialDictName_ = simpleMaterialChemistryModel::materialDict_.path() + "/" + simpleMaterialChemistryModel::materialDict_.name();
  const word pato_dir_ = getEnv("PATO_DIR");
  materialDictName_.replaceAll(pato_dir_, "$PATO_DIR");

  wordList validVATypes_;
  validVATypes_.append("FibrousMaterialTypeA");
  if (!foundInList(volumeAblationType_, validVATypes_)) {
    FatalErrorInFunction << "Only volume ablation type(s) valid in \"" <<  materialDictName_<< ".VolumeAblation\" : " << nl << validVATypes_ << exit(FatalError);
  }

  wordList validGPTypes_;
  validGPTypes_.append("FiniteRate");
  if (!foundInList(gasPropertiesType_, validGPTypes_)) {
    FatalErrorInFunction << "Only gas properties type(s) valid in \"" <<  materialDictName_<< ".GasProperties\" : " << nl << validGPTypes_ << exit(FatalError);
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    if (isA<VolumeAblationFvPatchScalarField>(T_.boundaryField()[patchI])) {
      volumeAblationPatches_.append(patchI);
    } else if (isA<speciesBCFvPatchScalarField>(Ydefault_.boundaryField()[patchI])) {
      volumeAblationPatches_.append(patchI);
    }
  }

  forAll(volumeAblationPatches_, i) {
    label patchI = volumeAblationPatches_[i];
    if (isA<VolumeAblationFvPatchScalarField>(T_.boundaryField()[patchI])) {
      VolumeAblationFvPatchScalarField& T_bf_ = refCast<VolumeAblationFvPatchScalarField>(T_.boundaryFieldRef()[patchI]);
      volumeAblationBC_.append(&T_bf_.VolumeAblationBoundaryConditions_);
    } else if (isA<speciesBCFvPatchScalarField>(Ydefault_.boundaryField()[patchI])) {
      speciesBCFvPatchScalarField& Ydefault_bf_ = refCast<speciesBCFvPatchScalarField>(Ydefault_.boundaryFieldRef()[patchI]);
      speciesBC_.append(&Ydefault_bf_.speciesBCBoundaryConditions_);
      Ydefault_.correctBoundaryConditions();
    }
  }

  // verify diffusion coefficients
  if(Dm_.size() != ns_mix) {
    FatalErrorInFunction << "Species diffusion coefficients size(" << Dm_.size() <<") != species mixture size(" << ns_mix << ")" << exit(FatalError);
  }

  // add the fields in the interpolation scheme table
  forAll(massFractions_, specI) {
    fields.add(Y_[specI]);
    fields.add(Dm_[specI]);
  }

  diffY_.resize(ns_mix);
  // Change write operator
  forAll(Y_, specI) {
    Y_[specI].writeOpt() = IOobject::NO_WRITE;
    word diffYName = "diffY["+Y_[specI].name()+"]";
    diffY_.set(
        specI,
        new volScalarField
        (
            IOobject
            (
                diffYName,
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Dm_[specI]*eps_g_*rho_g_/eta0_
        )
    );
    forAll(mesh_.boundaryMesh(), patchI) {
      if (isA<coupledMixedFvPatchScalarField>(Y_[specI].boundaryFieldRef()[patchI])) {
        coupledMixedFvPatchScalarField& cm_ =  refCast<coupledMixedFvPatchScalarField>(Y_[specI].boundaryFieldRef()[patchI]);
        cm_.setTnbrName(speciesNames_[specI]);
        cm_.setKappa(diffYName);

      }
    }
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- end --- Foam::speciesConservationMaterialChemistryModel::speciesConservationMaterialChemistryModel" << endl;
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SpeciesConservationMaterialChemistryModel::~SpeciesConservationMaterialChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SpeciesConservationMaterialChemistryModel::update()
{
  if(simpleMaterialChemistryModel::debug_) {
    Info<<"---  T_.correctBoundaryConditions() --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }
  T_.correctBoundaryConditions();
  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- updateYBC() --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }
  updateYBC();

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- epsphi, epsrho_g_, Tgi, pgi, rhogi --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }
  // Mass flux of pyrolysis gases
  tmp<surfaceScalarField> epsphi_tmp = mDotGFace_ & mesh_.Sf();
  const surfaceScalarField& epsphi = epsphi_tmp();

  // Volume-averaged gas density
  tmp<volScalarField> epsrho_g_tmp = eps_g_ * rho_g_;
  const volScalarField& epsrho_g_ = epsrho_g_tmp();

  tmp<scalarField> Tgi_tmp = T_.internalField();
  const scalarField& Tgi = Tgi_tmp();
  tmp<scalarField> pgi_tmp = p_.internalField();
  const scalarField& pgi = pgi_tmp();
  tmp<scalarField> rhogi_tmp = rho_g_.internalField();
  const scalarField& rhogi = rhogi_tmp();

  tmp<fv::convectionScheme<scalar> > mvConvection
  (
      fv::convectionScheme<scalar>::New
      (
          mesh_,
          fields,
          epsphi,
          mesh_.divScheme("div(epsphi,Yi)")
      )
  );

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- volumeAblationModel_.updateSolidCarbonMassFraction --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }
  volumeAblationModel_.updateSolidCarbonMassFraction(Y_[iCs_]);

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- chemistry_.solve --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }
  // computing rates and updating dtChem
  dtChem_ = chemistry_.solve
            (
                dtChem_,
                Tgi,
                pgi,
                rhogi
            );

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- omegaHeterogeneousRate_, omegaHeterogeneousEnergy_ --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }
  forAll(omegaHeterogeneousRate_, cellI) {
    omegaHeterogeneousRate_[cellI] = - chemistry_.RR(iCs_)[cellI]; // store heterogeneous reaction rate
  }
  omegaHeterogeneousRate_.max(0.0);

  // enthalpy of the solid lost with heterogeneous reactions
  omegaHeterogeneousEnergy_ = - h_c_ * omegaHeterogeneousRate_;

  Y_[iCs_] = 0.0 * Y_[iCs_]; // reset solid fraction to zero

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- solve Yi --- Foam::speciesConservationMaterialChemistryModel::update()" << endl;
  }


  tmp<volScalarField> Yt_tmp = 0.0 * Y_[0];
  volScalarField& Yt = const_cast<volScalarField&>(Yt_tmp());

  for (label i = 0; i < Y_.size(); i++) {
    if (i != iCs_) { //  Solid Carbon conservation is done in sMass.H (not transported)
      volScalarField& Yi = Y_[i];                  // Species i mass fraction field
      volScalarField& Dmi = Dm_[i];                // Species i average multicomponent diffusion coefficient field
      volScalarField& totalYpi = pyrolysisModel_.piTotal();     // Species i production rate by the pyrolysis reactions
      if (pi_.size()>0) {
        totalYpi = pi_[i];
      }

      diffY_[i] = Dmi / eta0_ * epsrho_g_;
      if (this->dynamicMesh_) {
        fvScalarMatrix YiEqn
        (
            fvm::ddt(epsrho_g_, Yi)                              // storage
            - fvm::div(fvc::interpolate(epsrho_g_)*mesh_.phi(), Yi) // mesh motion correction (ALE)
            + mvConvection->fvmDiv(epsphi, Yi)                   // convection
            - fvm::laplacian(Dm_[i]*eps_g_*rho_g_/eta0_, Yi)           // fickian diffusion
            - totalYpi                                           // species produced by the pyrolysis reactions
            ==
            eps_g_ * chemistry_.RR(i)			     // finite-rate chemistry production
        );

        YiEqn.solve("Yi");
      } else {
        fvScalarMatrix YiEqn
        (
            fvm::ddt(epsrho_g_, Yi)                              // storage
            + mvConvection->fvmDiv(epsphi, Yi)                   // convection
            - fvm::laplacian(Dm_[i]*eps_g_*rho_g_/eta0_, Yi)           // fickian diffusion
            - totalYpi                                           // species produced by the pyrolysis reactions
            ==
            eps_g_ * chemistry_.RR(i)			     // finite-rate chemistry production
        );

        YiEqn.solve("Yi");
      }

      Yi.max(0.0);
      Yt += Yi;
    }
  }

  for (label i = 0; i < Y_.size(); i++) {
    Y_[i] = Y_[i] / Yt;
  }



  forAll(massFractions_, specI) {
    forAll(massFractions_[specI], cellI) {
      massFractions_[specI][cellI] = Y_[specI][cellI];
      oldMassFractions_[specI][cellI] = Y_[specI].oldTime()[cellI];
    }
    forAll(mesh_.boundaryMesh(), patchI) {
      forAll(mesh_.boundaryMesh()[patchI], faceI) {
        massFractions_[specI].boundaryFieldRef()[patchI][faceI] = Y_[specI].boundaryField()[patchI][faceI];
        oldMassFractions_[specI].boundaryFieldRef()[patchI][faceI] = Y_[specI].oldTime().boundaryField()[patchI][faceI];
      }
    }
  }

  forAll(T_, cellI) {
    double * Z_Y_to_X_ = new double[ns_mix];
    forAll(massFractions_, specI) {
      Z_Y_to_X_[specI] = massFractions_[specI][cellI];
    }
    mix_().convert<Mutation::Thermodynamics::YE_TO_XE>(Z_Y_to_X_, Z_Y_to_X_);
    forAll(moleFractions_, specI) {
      moleFractions_[specI][cellI] = Z_Y_to_X_[specI];
    }
  }

  forAll(T_.boundaryField(), patchI) {
    forAll(T_.boundaryField()[patchI], faceI) {
      double * Z_Y_to_X_ = new double[ns_mix];

      forAll(massFractions_, specI) {
        Z_Y_to_X_[specI] = massFractions_[specI].boundaryField()[patchI][faceI];
      }
      mix_().convert<Mutation::Thermodynamics::YE_TO_XE>(Z_Y_to_X_, Z_Y_to_X_);
      forAll(moleFractions_, specI) {
        moleFractions_[specI].boundaryFieldRef()[patchI][faceI] = Z_Y_to_X_[specI];
      }
    }
  }
}

Switch Foam::SpeciesConservationMaterialChemistryModel::init()
{

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- init MaterialChemistry --- Foam::speciesConservationMaterialChemistryModel::init" << endl;
  }

  if (!isFile(constantPropertiesDictionary.path()+"/constantProperties")) {
    FatalErrorInFunction << "Unknown materialPropertiesDirectory " << constantPropertiesDictionary.path()
                         << exit(FatalError);
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- mixture --- Foam::speciesConservationMaterialChemistryModel::init" << endl;
  }

  optEq_.reset(new Mutation::MixtureOptions(mixtureMutation));
  optEq_().setStateModel("ChemNonEq1T");
  mix_.reset(new Mutation::Mixture(optEq_()));

  ns_mix=mix_().nSpecies();
  speciesNames_.resize(ns_mix);

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- speciesName --- Foam::speciesConservationMaterialChemistryModel::init" << endl;
  }

  if (speciesNames_.size()!=Y_.size()) {
    FatalErrorInFunction << "Mutation++ mixture has a different size than the CHEMKIN mixture." << exit(FatalError);
  }

  forAll(speciesNames_, specI) {
    speciesNames_[specI]=Y_[specI].name();
  }


  // Initialize the index from Mutation++
  speciesIndexMutation_.resize(speciesNames_.size());
  for (int i = 0; i < speciesNames_.size(); i++) {
    speciesIndexMutation_[i] = mix_->speciesIndex(speciesNames_[i]);
    if (speciesIndexMutation_[i]  < 0) {
      FatalErrorInFunction << "Mutation++ mixture does not match CHEMKIN mixture: " << nl << "Mutation++ species = ( ";
      for(int i = 0; i < mix_->nSpecies(); i++) {
        FatalErrorInFunction << word(mix_->speciesName(i)) << " ";
      }

      FatalErrorInFunction << ")" << nl <<"CHEMKIN species = ( ";
      forAll(Y_, i) {
        FatalErrorInFunction << Y_[i].name() << " ";
      }
      FatalErrorInFunction << ")" <<  exit(FatalError);
    }
  }

  moleFractions_.resize(ns_mix);
  massFractions_.resize(ns_mix);
  oldMassFractions_.resize(ns_mix);

  forAll(massFractions_, specI) {
    massFractions_.set
    (
        specI,
        new volScalarField
        (
            IOobject
            (
                "Y["+ speciesNames_[specI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Y_[specI]
        )
    );
    oldMassFractions_.set
    (
        specI,
        new volScalarField
        (
            IOobject
            (
                "Yold["+ speciesNames_[specI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Y_[specI]
        )
    );
    moleFractions_.set
    (
        specI,
        new volScalarField
        (
            IOobject
            (
                "X[" + speciesNames_[specI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Y_[specI]
        )
    );
  }

  // Extract indices of C(gr) for in-depth oxidation of the solid
  iCs_ = -1;
  for (int i = 0; i < Y_.size(); i++) {
    if (Y_[i].name() == "C(gr)") {
      iCs_ = i;
    }
  }

  if (iCs_<0) {
    FatalErrorInFunction << "C(gr) not found in the mixture." << exit(FatalError);
  }

  return true;
}

void Foam::SpeciesConservationMaterialChemistryModel::updateYBC()
{
  if (this->debug_) {
    Info << "--- start --- Foam::speciesConservationMaterialChemistrySolver<MaterialChemistryModel>::updateYBC()" << endl;
  }
  // Imposing fixedValues on the SEMB patches: Y_[j]
  forAll(volumeAblationPatches_, patchI) {
    scalarField rhoeUeCH_BF = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>("rhoeUeCH")).boundaryField()[patchI];
    scalarField blowingCorrection_BF = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>("blowingCorrection")).boundaryField()[patchI];
    const scalarField& invDx_BF    = mesh_.deltaCoeffs().boundaryField()[patchI];
    const scalarField& eps_g_BF     = eps_g_.boundaryField()[patchI];
    const scalarField& rho_g_BF     = rho_g_.boundaryField()[patchI];

    forAll(Y_, specI) {
      if (specI != iCs_) { //  Solid Carbon is a derived value
        scalarField& yBF = Y_[specI].boundaryFieldRef()[patchI];
        tmp<scalarField> Ysub_tmp = Y_[specI].boundaryField()[patchI].patchInternalField();
        const scalarField& Ysub = Ysub_tmp();
        const scalarField& Dmiw = Dm_[specI].boundaryField()[patchI];
        // update the gradient at the wall
        scalar Bim = 0;
        // updating BL values
        // note : here, we neglect advective flux in front of convective flux - jl, 7 Nov 2017
        forAll(rhoeUeCH_BF, faceI) {
          Bim = rhoeUeCH_BF[faceI]  * blowingCorrection_BF[faceI] / (Dmiw[faceI] * invDx_BF[faceI] * eps_g_BF[faceI] * rho_g_BF[faceI] / eta0_.value());
          if (Bim == 0) {
            FatalErrorInFunction << "Problem in the imposed value of the species mass fractions (updateYBC()). Bim = 0, rhoeUeCH_BF[faceI]  * blowingCorrection_BF[faceI] == 0" <<
                                 nl << "rhoeUeCH_BF=" << rhoeUeCH_BF << nl << "blowingCorrection_BF[faceI]" << blowingCorrection_BF[faceI] << exit(FatalError);
          }
          scalar Yie_ref=0;
          if (isA<VolumeAblationFvPatchScalarField>(T_.boundaryField()[patchI])) {
            Yie_ref = volumeAblationBC_[patchI]->Yie_ref[specI];
          } else if (isA<speciesBCFvPatchScalarField>(Ydefault_.boundaryField()[patchI])) {
            Yie_ref = speciesBC_[patchI]->Yie_ref[specI];
          }
          yBF[faceI] = 1 / (1 + 1/Bim) * Yie_ref + 1 / (1 + Bim) * Ysub[faceI];
        }
      }
    } // end loop on species
  }

  if (this->debug_) {
    Info << "--- end --- Foam::speciesConservationMaterialChemistrySolver<MaterialChemistryModel>::updateYBC()" << endl;
  }
}

// ************************************************************************* //
