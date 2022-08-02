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
    along with OpenFOAM.  If OnlyFiniteRatet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "OnlyFiniteRateMaterialChemistryModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OnlyFiniteRateMaterialChemistryModel::OnlyFiniteRateMaterialChemistryModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialChemistryModel(mesh, dictName),
mixtureMutation(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").lookup("mixture")),
mix_(simpleMaterialChemistryModel::mixture_),
elementNames_(simpleMaterialChemistryModel::elementNames_),
speciesNames_(simpleMaterialChemistryModel::speciesNames_),
massFractions_(simpleMaterialChemistryModel::massFractions_),
moleFractions_(simpleMaterialChemistryModel::moleFractions_),
T_(meshLookupOrConstructScalar(mesh,"Ta")),
p_(meshLookupOrConstructScalar(mesh,"p")),
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
pThermo_(psiReactionThermo::New(mesh)),
thermo_(pThermo_()),
pChemistry_(BasicFiniteRateChemistryModel<psiReactionThermo>::New(thermo_)),
chemistry_(pChemistry_()),
dtChem_(min(refCast<const BasicFiniteRateChemistryModel<psiReactionThermo>>(chemistry_).deltaTChem()[0],mesh.time().deltaT().value())),
composition_(thermo_.composition()),
Y_(composition_.Y()),
speciesIndexMutation_(simpleMaterialChemistryModel::speciesIndexMutation_),
init_(init()),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
rho_g_(gasPropertiesModel_.rho_g()),
eps_g_(gasPropertiesModel_.eps_g()),
Dm_(gasPropertiesModel_.Dm()),
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
eta0_(constantPropertiesDictionary.lookup("eta0"))
{
  if (Dm_.size()!=ns_mix) {
    FatalErrorInFunction << "Dm[i] not initialized." << exit(FatalError);
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
      if (isA<coupledMixedFvPatchScalarField>(Ydefault_.boundaryFieldRef()[patchI])) {
        coupledMixedFvPatchScalarField& cm_ =  refCast<coupledMixedFvPatchScalarField>(Y_[specI].boundaryFieldRef()[patchI]);
        cm_.setTnbrName(speciesNames_[specI]);
        cm_.setKappa(diffYName);
      }
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OnlyFiniteRateMaterialChemistryModel::~OnlyFiniteRateMaterialChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::OnlyFiniteRateMaterialChemistryModel::update()
{
  const tmp<scalarField> Tgi_tmp = T_.internalField();
  const scalarField& Tgi = Tgi_tmp();
  const tmp<scalarField> pgi_tmp = p_.internalField();
  const scalarField& pgi = pgi_tmp();
  const tmp<scalarField> rhogi_tmp = rho_g_.internalField();
  const scalarField& rhogi = rhogi_tmp();

  // computing rates and updating dtChem
  dtChem_ = chemistry_.solve
            (
                dtChem_,
                Tgi,
                pgi,
                rhogi
            );

  forAll(Y_, specI) {
    Y_[specI].correctBoundaryConditions();
  }

}

Switch Foam::OnlyFiniteRateMaterialChemistryModel::init()
{

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- init MaterialChemistry --- Foam::onlyFiniteRateMaterialChemistrySolver<MaterialChemistryModel>::onlyFiniteRateMaterialChemistrySolver" << endl;
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- mixture --- Foam::onlyFiniteRateMaterialChemistrySolver<MaterialChemistryModel>::onlyFiniteRateMaterialChemistrySolver" << endl;
  }

  optEq_.reset(new Mutation::MixtureOptions(mixtureMutation));
  optEq_().setStateModel("ChemNonEq1T");
  mix_.reset(new Mutation::Mixture(optEq_()));

  ns_mix=mix_().nSpecies();
  speciesNames_.resize(ns_mix);

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"--- speciesName --- Foam::onlyFiniteRateMaterialChemistrySolver<MaterialChemistryModel>::onlyFiniteRateMaterialChemistrySolver" << endl;
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

  forAll(massFractions_, specI) {
    massFractions_.set
    (
        specI,
        new volScalarField
        (
            IOobject
            (
                "Y[" + speciesNames_[specI]+"]",
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
                "X[" +speciesNames_[specI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            Y_[specI]
        )
    );
  }

  return true;
}


// ************************************************************************* //
