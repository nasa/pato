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
    const word& dictName,
    const label& currentPatchID,
    const dictionary dict
)
  :
mesh_(mesh),
phaseName_(phaseName),
dictName_(dictName),
currentPatchID_(currentPatchID),
dict_(dict),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
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
rhoeUeCH_(meshLookupOrConstructScalar(mesh, "rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0)))),
blowingCorrection_(meshLookupOrConstructScalar(mesh, "blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0)))),
materialDict_(
    IOobject
    (
        IOobject::groupName(dictName+"Properties", phaseName),
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
mixtureMutation(materialDict_.subDict("MaterialChemistry").lookup("mixture")),
materialPropertiesDirectory
(
    fileName(materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()
),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
rfFail(dimensionedScalar::lookupOrDefault("rfFail",constantPropertiesDictionary,dimensionSet(0, 1, 0, 0, 0, 0, 0),1e-6)),
debug_(materialDict_.lookupOrDefault<Switch>("debug","no"))
{
  optEq_ = new Mutation::MixtureOptions(mixtureMutation);
  optEq_->setStateModel("ChemNonEq1T");
  mix_ = new Mutation::Mixture(*optEq_);
  ns_mix = mix_->nSpecies();

  Yie_ref = new double [ns_mix];
  Info << "Reading boundary-layer edge species composition" << nl;
  for (int j = 0; j < ns_mix ; j++) {
    Yie_ref[j] =
        environmentDictionary_.lookupOrDefault<scalar>
        (
            "Yie[" + mix_->speciesName(j) + "]",
            0.0
        );

    Info << "Yie[" << mix_->speciesName(j) << "] = " <<  Yie_ref[j] << nl;
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::speciesBCBoundaryConditions::~speciesBCBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesBCBoundaryConditions::update()
{
  rhoeUeCH_.correctBoundaryConditions();
  blowingCorrection_.correctBoundaryConditions();
}



void Foam::speciesBCBoundaryConditions::write(Ostream& os) const
{
  os.writeKeyword("environmentDirectory") << environmentDirectory << token::END_STATEMENT << nl;
}


// ************************************************************************* //
