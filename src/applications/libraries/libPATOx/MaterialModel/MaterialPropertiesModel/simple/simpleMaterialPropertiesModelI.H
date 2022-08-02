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

#include "simpleMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::cp() const
{
  return const_cast<volScalarField&>(cp_);
}

inline Foam::volTensorField& Foam::simpleMaterialPropertiesModel::k() const
{
  return const_cast<volTensorField&>(k_);
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialPropertiesModel::solidRho()
{
  return solidRho_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialPropertiesModel::solidEps()
{
  return solidEps_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialPropertiesModel::solidEpsI() const
{
  return const_cast<PtrList<volScalarField>&>(solidEpsI_);
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialPropertiesModel::solidRhoI() const
{
  return const_cast<PtrList<volScalarField>&>(solidRhoI_);
}

inline Foam::volTensorField& Foam::simpleMaterialPropertiesModel::K() const
{
  return const_cast<volTensorField&>(K_);
}

inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::rho_s() const
{
  return const_cast<volScalarField&>(rho_s_);
}

inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::h_bar() const
{
  return const_cast<volScalarField&>(h_bar_);
}

inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::absorptivity() const
{
  return const_cast<volScalarField&>(absorptivity_);
}

inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::emissivity() const
{
  return const_cast<volScalarField&>(emissivity_);
}

inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::h_c() const
{
  return const_cast<volScalarField&>(h_c_);
}

inline Foam::volScalarField& Foam::simpleMaterialPropertiesModel::pyrolysisFlux() const
{
  return const_cast<volScalarField&>(pyrolysisFlux_);
}

inline Foam::scalar& Foam::simpleMaterialPropertiesModel::nSolidPhases() const
{
  return const_cast<scalar&>(nSolidPhases_);
}

// ************************************************************************* //

