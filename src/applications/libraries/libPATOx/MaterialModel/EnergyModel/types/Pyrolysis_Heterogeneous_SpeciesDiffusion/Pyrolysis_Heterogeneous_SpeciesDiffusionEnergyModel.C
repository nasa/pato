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
    along with OpenFOAM.  If Pyrolysis_Heterogeneous_SpeciesDiffusiont, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel::Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
PyrolysisEnergyModel(mesh, regionName),
materialChemistryModel_(refModel<simpleMaterialChemistryModel>()),
omegaHeterogeneousEnergy_(materialChemistryModel_.refVolField<scalar>("omegaHeterogeneousEnergy")),
Ediff_(gasPropertiesModel_.refVolField<scalar>("Ediff"))
{
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel::~Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel::update()
{

  PyrolysisEnergyModel::beforeSolve();

  if(this->dynamicMesh_) {
    // global energy balance
    solve
    (
        rho_s * cp * (fvm::ddt(T) - fvm::div(mesh_.phi(), T))         // storage - implicit in T, explicit in rhoCp - with mesh_ motion correction (ALE)
        + pyrolysisFlux_                                             // storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        + omegaHeterogeneousEnergy_                                   // storage (second part of the derivative - mass variation) - solid mass loss by heterogeneous reactions - explicit
        + fvc::ddt(epsgRhogEg) -  fvc::div(mesh_.phi(), epsgRhogEg)   // gas storage - explicit
        - fvm::laplacian(k, T)                                        // conduction - implicit in T, explicit in k
        - fvc::laplacian(GammaHg, p)                                  // convection - explicit
        - Ediff_                                                      // energy transported by diffusion of the species
    );
  } else {
    // global energy balance
    solve
    (
        rho_s * cp * (fvm::ddt(T))                                 // storage - implicit in T, explicit in rhoCp - with mesh motion correction (ALE)
        + pyrolysisFlux_                                           // storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        + omegaHeterogeneousEnergy_                                // storage (second part of the derivative - mass variation) - solid mass loss by heterogeneous reactions - explicit
        + fvc::ddt(epsgRhogEg)                                     // gas storage - explicit
        - fvm::laplacian(k, T)                                     // conduction - implicit in T, explicit in k
        - fvc::laplacian(GammaHg, p)                               // convection - explicit
        - Ediff_                                                   // energy transported by diffusion of the species
    );
  }

  PyrolysisEnergyModel::afterSolve();

}

// ************************************************************************* //
