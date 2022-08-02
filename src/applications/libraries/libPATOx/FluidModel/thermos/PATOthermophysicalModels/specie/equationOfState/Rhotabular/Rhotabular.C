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

#include "Rhotabular.H"
#include "IOstreams.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::Rhotabular<Specie>::Rhotabular(Istream& is)
  :
Specie(is)
{
  is.check("Rhotabular<Specie>::Rhotabular(Istream& is)");
  densityTable = interpolation2DTable<scalar>("constant/densityTable");
  densityTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);

}


template<class Specie>
Foam::Rhotabular<Specie>::Rhotabular(const dictionary& dict)
  :
Specie(dict),
densityTable(dict.subDict("equationOfState"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::Rhotabular<Specie>::write(Ostream& os) const
{
  return;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const Rhotabular<Specie>& ico)
{
  os  << static_cast<const Specie&>(ico);
  os.check("Ostream& operator<<(Ostream& os, const Rhotabular<Specie>& ico)");
  return os;
}


// ************************************************************************* //
