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

#include "simpleGasPropertiesModel.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


inline Foam::volScalarField& Foam::simpleGasPropertiesModel::M()
{
  return M_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::cp_g()
{
  return cp_g_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::mu()
{
  return mu_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::h_g()
{
  return h_g_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::rho_g()
{
  return rho_g_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::k_g()
{
  return k_g_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::eps_g()
{
  return eps_g_;
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleGasPropertiesModel::Dm()
{
  return Dm_;
}

inline Foam::volScalarField& Foam::simpleGasPropertiesModel::Ediff()
{
  return Ediff_;
}

// ************************************************************************* //

