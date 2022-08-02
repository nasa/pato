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
    along with OpenFOAM.  If ConstantEquilibriumt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ConstantEquilibriumMaterialChemistryModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ConstantEquilibriumMaterialChemistryModel::ConstantEquilibriumMaterialChemistryModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialChemistryModel(mesh, dictName),
mixtureMutation(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").lookup("mixture")),
mix_(simpleMaterialChemistryModel::mixture_),
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
moleFractions_(simpleMaterialChemistryModel::moleFractions_),
moleFractionGasInMaterial_(simpleMaterialChemistryModel::materialDict_.subDict("MaterialChemistry").template lookupOrDefault<List<Tuple2<word,scalar> > >("moleFractionGasInMaterial",List<Tuple2<word,scalar> >(0)))
{
  if (!isFile(constantPropertiesDictionary.path()+"/constantProperties")) {
    FatalErrorInFunction << "UnkConstantEquilibriumwn materialPropertiesDirectory " << constantPropertiesDictionary.path()
                         << exit(FatalError);
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"init MaterialChemistry" << endl;
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t mixture" << endl;
  }

  optEq_.reset(new Mutation::MixtureOptions(mixtureMutation));
  optEq_().setStateModel("Equil");
  mix_.reset(new Mutation::Mixture(optEq_()));

  ne_mix=mix_().nElements();
  elementNames_.resize(ne_mix);

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t elementName" << endl;
  }
  forAll(elementNames_, elemI) {
    elementNames_[elemI]=mix_().elementName(elemI);
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t Elements in mixture" << endl;
    Info << "ne = " << ne_mix << endl;
    forAll(elementNames_, elemI) {
      Info << elementNames_[elemI] << " ";
    }
    Info << endl;
  }

  moleFractions_.resize(ne_mix);
  massFractions_.resize(ne_mix);

  // Initialize Mutation++ to compute constantEquilibrium MaterialChemistry, thermodynamics and transport properties
  // Initialization of Mutation++ for constantEquilibrium MaterialChemistry
  p_Z.resize(ne_mix);      // mass fractions of elements
  p_Z0.resize(ne_mix);     // ref mass fractions of elements
  p_Zsub.resize(ne_mix);   // mass fractions of elements in the subsurface cell
  p_Zx.resize(ne_mix);     // mole fractions of elements

  if (moleFractionGasInMaterial_.size()>0) {
    if(simpleMaterialChemistryModel::debug_) {
      Info<<"\t moleFractionGasInMaterial_" << endl;
    }
    Info << "Use \"Zx[...]\" from \"moleFractionGasInMaterial\" instead of \"" << (word) constantPropertiesDictionary.path() << "/constantProperties\"." << endl;

    scalar sumZx_ = 0;
    scalarList values(moleFractionGasInMaterial_.size());
    wordList names(moleFractionGasInMaterial_.size());

    forAll(moleFractionGasInMaterial_, gasI) {
      names[gasI] = "Zx["+moleFractionGasInMaterial_[gasI].first()+"]";
      values[gasI] = moleFractionGasInMaterial_[gasI].second();
      sumZx_ += values[gasI];
    }
    if (sumZx_ != 1) {
      forAll(moleFractionGasInMaterial_, gasI) {
        values[gasI]/=sumZx_;
      }
    }
    forAll(moleFractionGasInMaterial_, gasI) {
      Info << names[gasI] << " = " << values[gasI] << endl;
      constantPropertiesDictionary.add(names[gasI] ,values[gasI]);
    }
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t read from constantPropertiesDictionary" << endl;
  }
  Info << "Reading initial elemental composition of the gas in the sample" << nl;
  scalar sumZx_ = 0;
  for (int i = 0; i < ne_mix; i++) {
    p_Zx[i] =
        constantPropertiesDictionary.lookupOrDefault<scalar>
        (
            "Zx[" + mix_().elementName(i) + "]",
            0.0
        );
    Info << "Zx[" << (word) mix_().elementName(i) << "] = " << p_Zx[i] << nl;
    p_Z[i] = p_Zx[i];
    sumZx_+=p_Zx[i];
  }

  if (sumZx_ != 1) {
    for (int i = 0; i < ne_mix; i++) {
      p_Zx[i]/=sumZx_;
      p_Z[i] = p_Zx[i];
    }
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t copyZ" << endl;
  }
  scalar copy_p_Z[ne_mix];
  forAll(p_Z, elemI) {
    copy_p_Z[elemI]=p_Z[elemI];
  }

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t mix XE_TO_YE" << endl;
  }
  mix_().convert<Mutation::Thermodynamics::XE_TO_YE>(copy_p_Z, copy_p_Z);  // convert in mass fractions

  if(simpleMaterialChemistryModel::debug_) {
    Info<<"\t massFractions_ and moleFractions_" << endl;
  }
  massFractions_.resize(elementNames_.size());
  moleFractions_.resize(elementNames_.size());
  forAll(massFractions_, elemI) {
    massFractions_.set
    (
        elemI,
        new volScalarField
        (
            IOobject
            (
                "Z["+ elementNames_[elemI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless, scalar(0))
        )
    );
    moleFractions_.set
    (
        elemI,
        new volScalarField
        (
            IOobject
            (
                "Zx["+ elementNames_[elemI]+"]",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless, scalar(0))
        )
    );
  }

  forAll(p_Z, elemI) {
    moleFractions_[elemI] == p_Zx[elemI];
    moleFractions_[elemI].boundaryFieldRef() ==  p_Zx[elemI];
    p_Z[elemI]=copy_p_Z[elemI];
    p_Z0[elemI] = p_Z[elemI];
    p_Zsub[elemI] = p_Z[elemI];
    massFractions_[elemI] == p_Z[elemI];
    massFractions_[elemI].boundaryFieldRef() ==  p_Z[elemI];
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ConstantEquilibriumMaterialChemistryModel::~ConstantEquilibriumMaterialChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ConstantEquilibriumMaterialChemistryModel::update()
{}

// ************************************************************************* //
