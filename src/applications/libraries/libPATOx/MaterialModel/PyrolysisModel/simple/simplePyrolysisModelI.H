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

#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


inline const Foam::fvMesh& Foam::simplePyrolysisModel::mesh() const
{
  return mesh_;
}

inline Foam::PtrList<volScalarField>& Foam::simplePyrolysisModel::piPyroReac()
{
  return const_cast<Foam::PtrList<volScalarField>&>(piPyroReac_);
}

inline Foam::PtrList<volScalarField>& Foam::simplePyrolysisModel::pi()
{
  return const_cast<Foam::PtrList<volScalarField>&>(pi_);
}

inline Foam::volScalarField& Foam::simplePyrolysisModel::piTotal()
{
  return const_cast<Foam::volScalarField&>(piTotal_);
}

inline Foam::volScalarField& Foam::simplePyrolysisModel::tau()
{
  return const_cast<Foam::volScalarField&>(tau_);
}

inline Foam::volScalarField& Foam::simplePyrolysisModel::Ta()
{
  return const_cast<Foam::volScalarField&>(Ta_);
}

inline Foam::volScalarField& Foam::simplePyrolysisModel::rho_s()
{
  return  const_cast<Foam::volScalarField&>(rho_s_);
}

inline Foam::volScalarField& Foam::simplePyrolysisModel::rho_v()
{
  return  const_cast<Foam::volScalarField&>(rho_v_);
}

inline Foam::volScalarField& Foam::simplePyrolysisModel::rho_c()
{
  return  const_cast<Foam::volScalarField&>(rho_c_);
}

inline Foam::PtrList<Foam::dimensionedScalar>& Foam::simplePyrolysisModel::hp()
{
  return  const_cast<Foam::PtrList<Foam::dimensionedScalar>&>(hp_);
}


// ************************************************************************* //

