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

    Class       : Foam::tabularTransport

    Description : Tabulated transport propertiy values for mu and kappa
		  read from tables.

    Note        :
                 - Dynamic viscosity evaluated to an expression in [Pa.s],
                   but internally uses [Pa.s/kmol].
                 - Thermal conductivity evaluated to an expression in [W/m/K],
                   but internally uses [W/m/K/kmol].

    SourceFiles :  tabularTransportI.H
                   tabularTransport.C

\*---------------------------------------------------------------------------*/

#ifndef tabularTransport_H
#define tabularTransport_H
#include "interpolation2DTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tabularTransport;

template<class Thermo>
inline tabularTransport<Thermo> operator+
(
    const tabularTransport<Thermo>&,
    const tabularTransport<Thermo>&
);

template<class Thermo>
inline tabularTransport<Thermo> operator-
(
    const tabularTransport<Thermo>&,
    const tabularTransport<Thermo>&
);

template<class Thermo>
inline tabularTransport<Thermo> operator*
(
    const scalar,
    const tabularTransport<Thermo>&
);

template<class Thermo>
inline tabularTransport<Thermo> operator==
(
    const tabularTransport<Thermo>&,
    const tabularTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabularTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class tabularTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabularTransport
:
public Thermo
{
  // Private data

  //- dynamic viscosity table reader and interpolator [kg/ms]
  interpolation2DTable<scalar> mu_;

  // thermal conductivity table reader and interpolator [W/(Km)]
  interpolation2DTable<scalar> kappa_;


  // Private Member Functions

  //- Construct from components
  inline tabularTransport
  (
      const Thermo& t,
      const interpolation2DTable<scalar> mu,
      const interpolation2DTable<scalar> kappa
  );


 public:

  // Constructors

  //- Construct as named copy
  inline tabularTransport(const word&, const tabularTransport&);

  //- Construct from Istream
  tabularTransport(Istream&);

  //- Construct from dictionary
  tabularTransport(const dictionary& dict);

  //- Construct and return a clone
  inline autoPtr<tabularTransport> clone() const;

  // Selector from Istream
  inline static autoPtr<tabularTransport> New(Istream& is);

  // Selector from dictionary
  inline static autoPtr<tabularTransport> New(const dictionary& dict);


  // Member functions

  //- Return the instantiated type name
  static word typeName() {
    return "tabular<" + Thermo::typeName() + '>';
  }

  //- Dynamic viscosity [kg/ms]
  inline scalar mu(const scalar T, const scalar p) const;

  //- Thermal conductivity [W/mK]
  inline scalar kappa(const scalar T, const scalar p) const;

  //- Thermal diffusivity of enthalpy [kg/ms]
  inline scalar alphah(const scalar T, const scalar p) const;

  // Species diffusivity (not implemented)
  //inline scalar D(const scalar p, const scalar T) const;

  //- Write to Ostream
  void write(Ostream& os) const;


  // Member operators

  inline tabularTransport& operator=(const tabularTransport&);
  inline void operator+=(const tabularTransport&);
  inline void operator-=(const tabularTransport&);
  inline void operator*=(const scalar);


  // Friend operators

  friend tabularTransport operator+ <Thermo>
  (
      const tabularTransport&,
      const tabularTransport&
  );

  friend tabularTransport operator- <Thermo>
  (
      const tabularTransport&,
      const tabularTransport&
  );

  friend tabularTransport operator* <Thermo>
  (
      const scalar,
      const tabularTransport&
  );

  friend tabularTransport operator== <Thermo>
  (
      const tabularTransport&,
      const tabularTransport&
  );


  // Ostream Operator

  friend Ostream& operator<< <Thermo>
  (
      Ostream&,
      const tabularTransport&
  );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabularTransportI.H"

#ifdef NoRepository
#   include "tabularTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
