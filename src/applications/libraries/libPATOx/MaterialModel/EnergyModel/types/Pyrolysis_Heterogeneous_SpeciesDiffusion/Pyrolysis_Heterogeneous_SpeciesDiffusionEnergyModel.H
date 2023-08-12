/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel

Description

SourceFiles
    Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel.H
    Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel_H
#define Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel_H

#include "PyrolysisEnergyModel.H"
#include "simpleMaterialChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel
:
public PyrolysisEnergyModel
{

 protected:

  //- Material Chemistry Model
  const simpleMaterialChemistryModel& materialChemistryModel_;

  //- Heterogeneous energy [W/m3] - updated in MaterialChemistryModel
  const volScalarField& omegaHeterogeneousEnergy_;

  //- Flux transported by diffusion of the species [W/m3]
  const volScalarField& Ediff_;

 public:

  //- Runtime type information
  TypeName("Pyrolysis_Heterogeneous_SpeciesDiffusion");

  // Constructors

  //- Construct from mesh and phase name
  Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~Pyrolysis_Heterogeneous_SpeciesDiffusionEnergyModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
