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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabularTransport<Thermo>::tabularTransport
(
    const Thermo& t,
    const interpolation2DTable<scalar> mu,
    const interpolation2DTable<scalar> kappa
)
  :
Thermo(t),
mu_(mu),
kappa_(kappa)
{
  mu_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
  kappa_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class Thermo>
inline Foam::tabularTransport<Thermo>::tabularTransport
(
    const word& name,
    const tabularTransport& ct
)
  :
Thermo(name, ct),
mu_(ct.mu_),
kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabularTransport<Thermo> >
Foam::tabularTransport<Thermo>::clone() const
{
  return autoPtr<tabularTransport<Thermo> >
         (
             new tabularTransport<Thermo>(*this)
         );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabularTransport<Thermo> >
Foam::tabularTransport<Thermo>::New
(
    Istream& is
)
{
  return autoPtr<tabularTransport<Thermo> >
         (
             new tabularTransport<Thermo>(is)
         );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabularTransport<Thermo> >
Foam::tabularTransport<Thermo>::New
(
    const dictionary& dict
)
{
  return autoPtr<tabularTransport<Thermo> >
         (
             new tabularTransport<Thermo>(dict)
         );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabularTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
  //Info << "mu value: " << mu_(T, p) << endl; // added to verify the well implementation of the table reader and interpolator;
  return mu_(T, p);
}


template<class Thermo>
inline Foam::scalar Foam::tabularTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
  // Info << "kappa value: " << kappa_(T, p) << endl; // added to verify the well implementation of the table reader and interpolator;
  return kappa_(T, p);
}


template<class Thermo>
inline Foam::scalar Foam::tabularTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
  return kappa(T, p) / this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabularTransport<Thermo>& Foam::tabularTransport<Thermo>::operator=
(
    const tabularTransport<Thermo>& pt
)
{
  Thermo::operator=(pt);

  mu_ = pt.mu_;
  kappa_ = pt.kappa_;

  return *this;

}

/* Operation removed due to compilation bug on Ubuntu 18, jl 22/02/2019 - please check usefulness
template<class Thermo>
inline void Foam::tabularTransport<Thermo>::operator+=
(
    const tabularTransport<Thermo>& pt
)
{
  scalar Y1 = this->Y();

  Thermo::operator+=(pt);

  if (mag(this->Y()) > SMALL) {
    Y1 /= this->Y();
    scalar Y2 = pt.Y()/this->Y();

    mu_ = Y1*mu_ + Y2*pt.mu_;
    kappa_ = Y1*kappa_ + Y2*pt.kappa_;
  }
}


template<class Thermo>
inline void Foam::tabularTransport<Thermo>::operator-=
(
    const tabularTransport<Thermo>& pt
)
{
  scalar Y1 = this->Y();

  Thermo::operator-=(pt);

  if (mag(this->Y()) > SMALL) {
    Y1 /= this->Y();
    scalar Y2 = pt.Y()/this->Y();

    mu_ = Y1*mu_ + Y2*pt.mu_;
    kappa_ = Y1*kappa_ + Y2*pt.kappa_;
  }
}
*/

template<class Thermo>
inline void Foam::tabularTransport<Thermo>::operator*=
(
    const scalar s
)
{
  Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabularTransport<Thermo> Foam::operator+
(
    const tabularTransport<Thermo>& pt1,
    const tabularTransport<Thermo>& pt2
)
{
  Thermo t
  (
      static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
  );

  if (mag(t.Y()) < SMALL) {
    return tabularTransport<Thermo>
           (
               t,
               0,
               pt1.mu_,
               pt1.kappa_
           );
  } else {
    scalar Y1 = pt1.Y()/t.Y();
    scalar Y2 = pt2.Y()/t.Y();

    return tabularTransport<Thermo>
           (
               t,
               Y1*pt1.mu_ + Y2*pt2.mu_,
               Y1*pt1.kappa_ + Y2*pt2.kappa_
           );
  }
}


template<class Thermo>
inline Foam::tabularTransport<Thermo> Foam::operator*
(
    const scalar s,
    const tabularTransport<Thermo>& pt
)
{
  return tabularTransport<Thermo>
         (
             s*static_cast<const Thermo&>(pt),
             pt.mu_,
             pt.kappa_
         );
}


// ************************************************************************* //
