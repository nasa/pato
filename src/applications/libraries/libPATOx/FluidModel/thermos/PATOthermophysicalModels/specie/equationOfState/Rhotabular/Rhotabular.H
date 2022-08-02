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

    Class       :    Foam::Rhotabular

    Description :    Tabulated form of equation of state, reading the density
                     values from tables "densityTable".

    SourceFiles :    RhotabularI.H
                     Rhotabular.C

\*---------------------------------------------------------------------------*/

#ifndef Rhotabular_H
#define Rhotabular_H

#include "autoPtr.H"
#include "interpolation2DTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

namespace Foam
{

template<class Specie>
class Rhotabular;

template<class Specie>
Rhotabular<Specie> operator+
(
    const Rhotabular<Specie>&,
    const Rhotabular<Specie>&
);

template<class Specie>
Rhotabular<Specie> operator-
(
    const Rhotabular<Specie>&,
    const Rhotabular<Specie>&
);

template<class Specie>
Rhotabular<Specie> operator*
(
    const scalar,
    const Rhotabular<Specie>&
);

template<class Specie>
Rhotabular<Specie> operator==
(
    const Rhotabular<Specie>&,
    const Rhotabular<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const Rhotabular<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class Rhotabular Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Rhotabular
:
public Specie
{
  // Private data

  //- Density table reader and interpolator :
  //  Note: input in [kg/m3], but internally uses [kg/m3/kmol]
  interpolation2DTable<scalar> densityTable;


 public:

  // Constructors

  //- Construct from components
  inline Rhotabular
  (
      const Specie& sp,
      const interpolation2DTable<scalar>& densityTable
  );

  //- Construct from Istream
  Rhotabular(Istream&);

  //- Construct from dictionary
  Rhotabular(const dictionary& dict);

  //- Construct as named copy
  inline Rhotabular(const word& name, const Rhotabular&);

  //- Construct and return a clone
  inline autoPtr<Rhotabular> clone() const;

  // Selector from Istream
  inline static autoPtr<Rhotabular> New(const dictionary& dict);


  // Member functions

  //- Return the instantiated type name
  static word typeName() {
    return "Rhotabular<" + word(Specie::typeName_()) + '>';
  }


  // Fundamental properties

  //- Is the equation of state is incompressible i.e. rho != f(p)
  static const bool incompressible = false;

  //- Is the equation of state is isochoric i.e. rho = const
  static const bool isochoric = false;

  //- Return density [kg/m^3]
  inline scalar rho(scalar p, scalar T) const;

  //- Return enthalpy departure [J/kg]
  inline scalar H(scalar p, scalar T) const;

  //- Return sensible energy [J/kg]
  inline scalar Es(scalar p, scalar T) const;

  //- Return Cp departure [J/(kg K]
  inline scalar Cp(scalar p, scalar T) const;

  //- Return Cv departure [J/(kg K]
  inline scalar Cv(scalar p, scalar T) const;

  //- Return entropy [J/(kg K)]
  inline scalar S(const scalar p, const scalar T) const;

  //- Return compressibility rho/p [s^2/m^2]
  inline scalar psi(scalar p, scalar T) const;

  //- Return compression factor []
  inline scalar Z(scalar p, scalar T) const;

  //- Return (Cp - Cv) [J/(kg K]
  inline scalar CpMCv(scalar p, scalar T) const;


  // IO

  //- Write to Ostream
  void write(Ostream& os) const;


  // Member operators

  inline void operator+=(const Rhotabular&);
  inline void operator-=(const Rhotabular&);
  inline void operator*=(const scalar);


  // Friend operators

  friend Rhotabular operator+ <Specie>
  (
      const Rhotabular&,
      const Rhotabular&
  );

  friend Rhotabular operator- <Specie>
  (
      const Rhotabular&,
      const Rhotabular&
  );

  friend Rhotabular operator* <Specie>
  (
      const scalar s,
      const Rhotabular&
  );

  friend Rhotabular operator== <Specie>
  (
      const Rhotabular&,
      const Rhotabular&
  );


  // Ostream Operator

  friend Ostream& operator<< <Specie>
  (
      Ostream&,
      const Rhotabular&
  );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RhotabularI.H"

#ifdef NoRepository
#   include "Rhotabular.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
