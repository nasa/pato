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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const EquationOfState& pt,
    const scalar Hf,
    const interpolation2DTable<scalar> cpTable,
    const interpolation2DTable<scalar> hTable
)
  :
EquationOfState(pt),
Hf_(Hf),
cpTable(cpTable),
hTable(hTable)

{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const hTabularThermo& pt
)
  :
EquationOfState(pt),
Hf_(pt.Hf_),
cpTable(pt.cpTable),
hTable(pt.hTable)
{}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const word& name,
    const hTabularThermo& pt
)
  :
EquationOfState(name, pt),
Hf_(pt.Hf_),
cpTable(pt.cpTable),
hTable(pt.hTable)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::limit
(
    const scalar T
) const
{
  return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::Cp
(
    const scalar p, const scalar T
) const
{
  //Info << "cp for T= " << T << ", p= " << p  << ": " << cpTable(T, p) << endl;
  return cpTable(p, T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
  return hTable(p, T) + EquationOfState::H(p, T);

}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
  return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::Hc()
const
{
  return Hf_;
}

template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
  notImplemented
  (
      "scalar hTabularThermo<EquationOfState>::"
      "s(const scalar p, const scalar T) const"
  );
  return T;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>&
Foam::hTabularThermo<EquationOfState>::operator=
(
    const hTabularThermo<EquationOfState>& pt
)
{
  EquationOfState::operator=(pt);

  Hf_ = pt.Hf_;
  cpTable = pt.cpTable;
  hTable = pt.hTable;

  return *this;
}

/* Operation removed due to compilation bug on Ubuntu 18, jl 22/02/2019 - please check usefulness
template<class EquationOfState>
inline void Foam::hTabularThermo<EquationOfState>::operator+=
(
    const hTabularThermo<EquationOfState>& pt
)
{
  scalar W1 = this->W();

  EquationOfState::operator+=(pt);

  if (mag(this->W()) > SMALL) {
    W1 /= this->W();
    const scalar W2 = pt.W()/this->W();

    Hf_ = W1*Hf_ + W2*pt.Hf_;
    cpTable = W1*cpTable + W2*pt.cpTable;
    hTable = W1*hTable + W2*pt.hTable;
  }
}
*/

template<class EquationOfState>
inline void Foam::hTabularThermo<EquationOfState>::operator*=
(
    const scalar s
)
{
  EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator+
(
    const hTabularThermo<EquationOfState>& pt1,
    const hTabularThermo<EquationOfState>& pt2
)
{
  EquationOfState eofs = pt1;
  eofs += pt2;

  if (mag(eofs.W()) < SMALL) {
    return hTabularThermo<EquationOfState>
           (
               eofs,
               pt1.Hf_,
               pt1.cpTable,
               pt1.hTable
           );
  }
  {
    const scalar W1 = pt1.W()/eofs.W();
    const scalar W2 = pt2.W()/eofs.W();

    return hTabularThermo<EquationOfState>
           (
               eofs,
               W1*pt1.Hf_ + W2*pt2.Hf_,
               W1*pt1.cpTable + W2*pt2.cpTable,
               W1*pt1.hTable + W2*pt2.hTable
           );
  }
}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hTabularThermo<EquationOfState>& pt
)
{
  return hTabularThermo<EquationOfState>
         (
             s*static_cast<const EquationOfState&>(pt),
             pt.Hf_,
             pt.cpTable,
             pt.hTable
         );
}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator==
(
    const hTabularThermo<EquationOfState>& pt1,
    const hTabularThermo<EquationOfState>& pt2
)
{
  EquationOfState eofs
  (
      static_cast<const EquationOfState&>(pt1)
      == static_cast<const EquationOfState&>(pt2)
  );

  const scalar W1 = pt1.W()/eofs.W();
  const scalar W2 = pt2.W()/eofs.W();

  return hTabularThermo<EquationOfState>
         (
             eofs,
             W2*pt2.Hf_     - W1*pt1.Hf_,
             W2*pt2.cpTable - W1*pt1.cpTable,
             W2*pt2.hTable  - W1*pt1.hTable
         );
}


// ************************************************************************* //
