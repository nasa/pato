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

#include "DarcyLaw_HeterogeneousMassModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DarcyLaw_HeterogeneousMassModel::DarcyLaw_HeterogeneousMassModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
DarcyLawMassModel(mesh,dictName),
materialChemistryModel_(refModel<simpleMaterialChemistryModel>()),
omegaHeterogeneousRate_(materialChemistryModel_.refVolField<scalar>("omegaHeterogeneousRate"))
{
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DarcyLaw_HeterogeneousMassModel::~DarcyLaw_HeterogeneousMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DarcyLaw_HeterogeneousMassModel::update()
{
  // semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
  // momentum conservation- substituted inside the mass conservation equation)

  DarcyLawMassModel::beforeSolve();

  // compute implicitly the new pressure
  if(!simpleMassModel::dynamicMesh_) {
    solve
    (
        fvm::ddt(Eta, p)                                  // storage
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
        - omegaHeterogeneousRate_                         // source : heterogeneous reactions
    );
  } else {
    solve
    (
        fvm::ddt(Eta, p)                                  // storage
        - fvm::div(fvc::interpolate(Eta)*mesh_.phi(), p)   // mesh motion correction (ALE)
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
        - omegaHeterogeneousRate_                          // source : heterogeneous reactions
    );
  }

  DarcyLawMassModel::afterSolve();
}


// ************************************************************************* //
