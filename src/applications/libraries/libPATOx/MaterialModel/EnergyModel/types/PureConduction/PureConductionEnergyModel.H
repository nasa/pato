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
    along with OpenFOAM.  If PureConductiont, see <http://www.gnu.org/licenses/>.

Class
    Foam::PureConductionEnergyModel

Description

SourceFiles
    PureConductionEnergyModel.H
    PureConductionEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef PureConductionEnergyModel_H
#define PureConductionEnergyModel_H

#include "simpleEnergyModel.H"
#include "simpleMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PureConductionEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class PureConductionEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Mesh reference
  const fvMesh& mesh_;

  //- Dictionary name
  const word& dictName_;

  //- Material properties model
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- Temperature [K]
  volScalarField& T_;

  //- Solid mass density [kg/m3]
  volScalarField& rho_s_;

  //- Solid heat capacity [J/kg/K]
  volScalarField& cp_;

  //- Solid thermal conductivity [W/m/K]
  volTensorField& k_;

 public:

  //- Runtime type information
  TypeName("PureConduction");

  // Constructors

  //- Construct from mesh and phase name
  PureConductionEnergyModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PureConductionEnergyModel();

  //- Member Functions

  //- Update the energy conservation fields
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef PureConductionRepository
#include "PureConductionEnergyModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
