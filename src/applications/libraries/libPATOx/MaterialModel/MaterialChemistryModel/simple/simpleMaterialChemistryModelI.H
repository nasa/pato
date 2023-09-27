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


inline const autoPtr<Mutation::Mixture>& Foam::simpleMaterialChemistryModel::mixture() const
{
  return mixture_;
}

inline const Foam::labelList& Foam::simpleMaterialChemistryModel::speciesIndexMutation() const
{
  return  const_cast<labelList&>(speciesIndexMutation_);
}

inline const Foam::wordList& Foam::simpleMaterialChemistryModel::elementNames() const
{
  return const_cast<wordList&>(elementNames_);
}

inline const Foam::wordList& Foam::simpleMaterialChemistryModel::speciesNames() const
{
  return const_cast<wordList&>(speciesNames_);
}

inline const Foam::PList<Foam::volScalarField>& Foam::simpleMaterialChemistryModel::massFractions() const
{
  return const_cast<PList<volScalarField>&>(massFractions_);
}

inline const Foam::PList<Foam::volScalarField>& Foam::simpleMaterialChemistryModel::oldMassFractions() const
{
  return const_cast<PList<volScalarField>&>(oldMassFractions_);
}

inline const Foam::PList<Foam::volScalarField>& Foam::simpleMaterialChemistryModel::moleFractions() const
{
  return const_cast<PList<volScalarField>&>(moleFractions_);
}

inline const Foam::scalar& Foam::simpleMaterialChemistryModel::dtChem() const
{
  return const_cast<scalar&>(dtChem_);
}

// ************************************************************************* //

