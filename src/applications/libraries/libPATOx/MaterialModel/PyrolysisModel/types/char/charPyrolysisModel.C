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
    along with OpenFOAM.  If LinearArrheniust, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "charPyrolysisModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::charPyrolysisModel::charPyrolysisModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
virginPyrolysisModel(mesh, regionName),
startModelInit_(startModelInit())
{
  if (this->debug_) {
    Info << getTabLevel() << "debug: start --- Foam::charPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  tau*=0; // char state
  // Solid densities update
  forAll(solidRho_, phaseI) {
    solidRho_[phaseI] = solidRhoI_[phaseI] * oneVolScalarField ;
    for(int i = 0; i < sizeXsi; i++) {
      if (hashXsi[i][0]-1 == phaseI) {
        solidRho_[phaseI] -= solidRhoI_[phaseI] * Fp[i];
      }
    }
  }

  // this model is initialized
  modelInitialized();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::charPyrolysisModel::~charPyrolysisModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::charPyrolysisModel::update()
{}

// ************************************************************************* //
