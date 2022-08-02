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

#include "simpleMaterialChemistryModel.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


inline autoPtr<Mutation::Mixture>& Foam::simpleMaterialChemistryModel::mixture()
{
  return mixture_;
}

inline Foam::labelList& Foam::simpleMaterialChemistryModel::speciesIndexMutation() const
{
  return  const_cast<labelList&>(speciesIndexMutation_);
}

inline Foam::wordList& Foam::simpleMaterialChemistryModel::elementNames() const
{
  return const_cast<wordList&>(elementNames_);
}

inline Foam::wordList& Foam::simpleMaterialChemistryModel::speciesNames() const
{
  return const_cast<wordList&>(speciesNames_);
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialChemistryModel::massFractions() const
{
  return const_cast<PtrList<volScalarField>&>(massFractions_);
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialChemistryModel::oldMassFractions() const
{
  return const_cast<PtrList<volScalarField>&>(oldMassFractions_);
}

inline Foam::PtrList<Foam::volScalarField>& Foam::simpleMaterialChemistryModel::moleFractions() const
{
  return const_cast<PtrList<volScalarField>&>(moleFractions_);
}

inline Foam::scalar& Foam::simpleMaterialChemistryModel::dtChem() const
{
  return const_cast<scalar&>(dtChem_);
}

// ************************************************************************* //

