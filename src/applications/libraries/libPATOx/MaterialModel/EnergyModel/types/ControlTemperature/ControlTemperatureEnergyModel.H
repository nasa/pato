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
    Foam::ControlTemperatureEnergyModel

Description
    Fixed temperature in the domain function of time
    T(t) = coeff1 + coeff2 * t
    The coefficients a and b are given as input in the
    materialproperties file. Default values are 300 and 0

SourceFiles
    ControlTemperatureEnergyModel.H
    ControlTemperatureEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef ControlTemperatureEnergyModel_H
#define ControlTemperatureEnergyModel_H

#include "simpleEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ControlTemperatureEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class ControlTemperatureEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Temperature [K]
  volScalarField& T_;

  //- Solid mass density [kg/m3]
  volScalarField& rho_s_;

  //- T-t function coefficents
  dimensionedScalar coeff1;
  dimensionedScalar coeff2;

 public:

  //- Runtime type information
  TypeName("ControlTemperature");

  // Constructors

  //- Construct from mesh and phase name
  ControlTemperatureEnergyModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~ControlTemperatureEnergyModel();

  //- Member Functions

  //- Update the energy conservation fields
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef ControlTemperatureRepository
#include "ControlTemperatureEnergyModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
