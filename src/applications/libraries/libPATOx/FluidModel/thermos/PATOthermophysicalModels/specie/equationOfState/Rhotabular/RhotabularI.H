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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::Rhotabular<Specie>::Rhotabular
(
    const Specie& sp,
    const interpolation2DTable<scalar>& densityTable
)
  :
Specie(sp),
densityTable(densityTable)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Rhotabular<Specie>::Rhotabular
(
    const word& name,
    const Rhotabular<Specie>& ico
)
  :
Specie(name, ico),
densityTable(ico.densityTable)
{}


template<class Specie>
inline Foam::autoPtr<Foam::Rhotabular<Specie> >
Foam::Rhotabular<Specie>::clone() const
{
  return autoPtr<Rhotabular<Specie>>
         (
             new Rhotabular<Specie>(*this)
         );
}


template<class Specie>
inline Foam::autoPtr<Foam::Rhotabular<Specie> >
Foam::Rhotabular<Specie>::New(const dictionary& dict)
{
  return autoPtr<Rhotabular<Specie>>
         (
             new Rhotabular<Specie>(dict)
         );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::rho
(
    scalar p,
    scalar T
) const
{
  return densityTable(T, p);
}

// Could be deleted while returning zero, but we kept
// the same structure the icoPolynomialI.H file is defined.

template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::Es
(
    scalar p,
    scalar T
) const
{
  return 0;
}


template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::H
(
    scalar p,
    scalar T
) const
{
  return 0;
}

template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
  return 0;
}

template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::Cv
(
    scalar p,
    scalar T
) const
{
  return 0;
}

template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::S
(
    scalar p,
    scalar T
) const
{
  return 0;
}


template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::psi
(
    scalar p,
    scalar T
) const
{

  return (rho(p, T) / (p + SMALL));
}


template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::Z
(
    scalar p,
    scalar T
) const
{
  return 0;
}

//Not Sure , need to be verified :
//--------------------------------------------------------------
template<class Specie>
inline Foam::scalar Foam::Rhotabular<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
  return -(p/sqr(densityTable(T,p)))*(densityTable(T,p)/T);
}
//--------------------------------------------------------------

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::Rhotabular<Specie>::operator+=(const Rhotabular<Specie>& ico)
{
  Specie::operator+=(ico);

}


template<class Specie>
inline void Foam::Rhotabular<Specie>::operator-=(const Rhotabular<Specie>& ico)
{
  Specie::operator-=(ico);

}


template<class Specie>
inline void Foam::Rhotabular<Specie>::operator*=(const scalar s)
{
  Specie::operator*=(s);

}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
Foam::Rhotabular<Specie> Foam::operator+
(
    const Rhotabular<Specie>& ico1,
    const Rhotabular<Specie>& ico2
)
{
  Specie sp
  (
      static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2)
  );

  if (mag(sp.Y()) < SMALL) {
    return Rhotabular<Specie>
           (
               sp,
               ico1.densityTable
           );
  } else {
    const scalar Y1 = ico1.Y()/sp.Y();
    const scalar Y2 = ico2.Y()/sp.Y();

    return Rhotabular<Specie>
           (
               sp,
               Y1*ico1.densityTable + Y2*ico2.densityTable
           );
  }
}


template<class Specie>
Foam::Rhotabular<Specie> Foam::operator-
(
    const Rhotabular<Specie>& ico1,
    const Rhotabular<Specie>& ico2
)
{
  Specie sp
  (
      static_cast<const Specie&>(ico1)
      - static_cast<const Specie&>(ico2)
  )
  ;

  if (mag(sp.Y()) < SMALL) {
    return Rhotabular<Specie>
           (
               sp,
               ico1.densityTable
           );
  } else {
    const scalar Y1 = ico1.Y()/sp.Y();
    const scalar Y2 = ico2.Y()/sp.Y();

    return Rhotabular<Specie>
           (
               sp,
               Y1*ico1.densityTable - Y2*ico2.densityTable
           );
  }

}


template<class Specie>
Foam::Rhotabular<Specie> Foam::operator*
(
    const scalar s,
    const Rhotabular<Specie>& ico
)
{
  return Rhotabular<Specie>
         (
             s*static_cast<const Specie&>(ico),
             ico.densityTable
         );
}


template<class Specie>
Foam::Rhotabular<Specie> Foam::operator==
(
    const Rhotabular<Specie>& ico1,
    const Rhotabular<Specie>& ico2
)
{
  Specie sp
  (
      static_cast<const Specie&>(ico1)
      == static_cast<const Specie&>(ico2)
  );

  const scalar Y1 = ico1.Y()/sp.Y();
  const scalar Y2 = ico2.Y()/sp.Y();

  return Rhotabular<Specie>
         (
             sp,
             Y2*ico2.densityTable - Y1*ico1.densityTable
         );
}


// ************************************************************************* //
