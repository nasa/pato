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

#include "tabularTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::tabularTransport<Thermo>::tabularTransport(Istream& is)
  :
Thermo(is)
{
  is.check("tabularTransport::tabularTransport(Istream& is)");
  mu_ = interpolation2DTable<scalar>("constant/muTable");
  mu_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
  kappa_ = interpolation2DTable<scalar>("constant/kappaTable");
  kappa_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class Thermo>
Foam::tabularTransport<Thermo>::tabularTransport(const dictionary& dict)
  :
Thermo(dict),
mu_(dict.subDict("transport").subDict("mu")),
kappa_(dict.subDict("transport").subDict("kappa"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::tabularTransport<Thermo>::tabularTransport::write(Ostream& os) const
{
  os  << this->name() << endl;
  os  << token::BEGIN_BLOCK  << incrIndent << nl;

  Thermo::write(os);

  dictionary dict("transport");
  os  << indent << dict.dictName() << dict;

  os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const tabularTransport<Thermo>& ct)
{
  operator<<(os, static_cast<const Thermo&>(ct));

  os.check("Ostream& operator<<(Ostream&, const tabularTransport&)");

  return os;
}


// ************************************************************************* //
