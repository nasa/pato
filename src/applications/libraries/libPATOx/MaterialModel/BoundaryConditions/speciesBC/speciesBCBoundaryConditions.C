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
    along with OpenFOAM.  If Bprimet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "speciesBCBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesBCBoundaryConditions::speciesBCBoundaryConditions
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& regionName,
    const label& currentPatchID,
    const dictionary dict
)
  :
mesh_(mesh),
phaseName_(phaseName),
regionName_(regionName),
currentPatchID_(currentPatchID),
dict_(dict),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName)),
environmentDirectory(fileName(dict_.lookup("environmentDirectory")).expand()),
environmentDictionary_
(
    IOobject
    (
        "environmentComposition",
        environmentDirectory,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mixtureMutation(energyModel_.materialDict().subDict("MaterialChemistry").lookup("mixture")),
rfFail(energyModel_.createDimScalarProp("rfFail",true,dimensionedScalar("0",dimensionSet(0, 1, 0, 0, 0, 0, 0),1e-6))),
debug_(energyModel_.materialDict().lookupOrDefault<Switch>("debug","no"))
{
  // Create new fields in Energy Model
  scalarFields_.insert("Ta",energyModel_.createVolFieldIfNotFound<scalar>(energyModel_,"Ta"));
  scalarFields_.insert("rho_s",energyModel_.createVolFieldIfNotFound<scalar>(energyModel_,"rho_s"));
  scalarFields_.insert("rhoeUeCH",energyModel_.createVolFieldIfNotFound<scalar>(energyModel_,"rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0))));
  scalarFields_.insert("blowingCorrection",energyModel_.createVolFieldIfNotFound<scalar>(energyModel_,"blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0))));

  optEq_ = new Mutation::MixtureOptions(mixtureMutation);
  optEq_->setStateModel("ChemNonEq1T");
  mix_ = new Mutation::Mixture(*optEq_);
  ns_mix = mix_->nSpecies();

  Yie_ref = new double [ns_mix];
  Info << simpleModel::getTabLevel() << "Reading boundary-layer edge species composition" << nl;
  for (int j = 0; j < ns_mix ; j++) {
    Yie_ref[j] =
        environmentDictionary_.lookupOrDefault<scalar>
        (
            "Yie[" + mix_->speciesName(j) + "]",
            0.0
        );

    Info << simpleModel::getTabLevel() << "Yie[" << mix_->speciesName(j) << "] = " <<  Yie_ref[j] << nl;
  }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::speciesBCBoundaryConditions::~speciesBCBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesBCBoundaryConditions::update()
{
  scalarFields_["rhoeUeCH"].correctBoundaryConditions();
  scalarFields_["blowingCorrection"].correctBoundaryConditions();
}



void Foam::speciesBCBoundaryConditions::write(Ostream& os) const
{
  fileName envDir = environmentDirectory;
  os.writeKeyword("environmentDirectory") << envDir.replaceAll(getEnv("PATO_DIR"),"$PATO_DIR") << token::END_STATEMENT << nl;
}


// ************************************************************************* //
