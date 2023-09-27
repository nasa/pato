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
    along with OpenFOAM.  If ConstantFiniteRatet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ConstantFiniteRateMaterialChemistryModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ConstantFiniteRateMaterialChemistryModel::ConstantFiniteRateMaterialChemistryModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleMaterialChemistryModel(mesh, regionName),
mixtureMutationFiniteRate(materialDict_.subDict("MaterialChemistry").lookup("mixtureMutationFiniteRate")),
mixFiniteRate(simpleMaterialChemistryModel::mixture_),
elementNames_(simpleMaterialChemistryModel::elementNames_),
speciesNames_(simpleMaterialChemistryModel::speciesNames_),
massFractions_(simpleMaterialChemistryModel::massFractions_),
moleFractions_(simpleMaterialChemistryModel::moleFractions_)
{
  optsFiniteRate.reset(new Mutation::MixtureOptions(mixtureMutationFiniteRate));
  optsFiniteRate().setStateModel("ChemNonEq1T");
  mixFiniteRate.reset(new Mutation::Mixture(optsFiniteRate()));
  ns_mix = mixFiniteRate().nSpecies();
  speciesNames_.resize(ns_mix);
  forAll(speciesNames_, specieI) {
    speciesNames_[specieI]=mixFiniteRate().speciesName(specieI);
  }

  massFractions_.resize(speciesNames_.size());
  moleFractions_.resize(speciesNames_.size());

  forAll(massFractions_, specieI) {
    massFractions_.set(specieI,createVolField<scalar>("Y["+ speciesNames_[specieI]+"]",dimensionedScalar("0", dimless, scalar(0))));
    moleFractions_.set(specieI,createVolField<scalar>("X["+ speciesNames_[specieI]+"]",dimensionedScalar("0", dimless, scalar(0))));
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ConstantFiniteRateMaterialChemistryModel::~ConstantFiniteRateMaterialChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ConstantFiniteRateMaterialChemistryModel::update()
{}

// ************************************************************************* //
