/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 PATO-community
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Based on tilasoldo, Yuusha and chriss85 contribution to OpenFOAM, this new
    thermophysical model has been modified and checked by PATO-community.

    The Interpolation function used in this updated file is that of OpenFoam
    called "interpolation2DTable".

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

#include "hTabularThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    Istream& is
)
  :
EquationOfState(is),
Hf_(readScalar(is))
{
  Hf_ *= this->W();
  cpTable = interpolation2DTable<scalar>("constant/cpTable");
  hTable = interpolation2DTable<scalar>("constant/hTable");
  cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
  hTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class EquationOfState>
Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const dictionary& dict
)
  :
EquationOfState(dict),
Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf"))),
cpTable(dict.subDict("thermodynamics").subDict("Cp")),
hTable(dict.subDict("thermodynamics").subDict("h"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabularThermo<EquationOfState>::write
(
    Ostream& os
) const
{
  EquationOfState::write(os);

  dictionary dict("thermodynamics");
  dict.add("Hf", Hf_);
  os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hTabularThermo<EquationOfState>& pt
)
{
  os  << static_cast<const EquationOfState&>(pt) << tab
      << pt.Hf_ << tab;

  os.check
  (
      "operator<<"
      "("
      "Ostream&, "
      "const hTabularThermo<EquationOfState>&"
      ")"
  );

  return os;
}


// ************************************************************************* //
