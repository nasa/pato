/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

Typedef
    Foam::compressible::turbulenceModel

Typedef
    Foam::compressible::RASModel

Typedef
    Foam::compressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for compressible flow
    based on the standard laminar transport package.

SourceFiles
    turbulentUserThermo.C
    turbulentUserThermos.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentUserThermo_H
#define turbulentUserThermo_H

#include "CompressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "basicUserThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressibleUser
{
typedef ThermalDiffusivity<CompressibleTurbulenceModel<basicUserThermo>>
    turbulenceModel;

typedef laminarModel<turbulenceModel> laminarModel;
typedef RASModel<EddyDiffusivity<turbulenceModel>> RASModel;
typedef LESModel<EddyDiffusivity<turbulenceModel>> LESModel;

template<class BasicCompressibleTurbulenceModel>
autoPtr<BasicCompressibleTurbulenceModel> New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const typename BasicCompressibleTurbulenceModel::transportModel&
    transport,
    const word& propertiesName = turbulenceModel::propertiesName
);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "turbulentUserThermo.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
