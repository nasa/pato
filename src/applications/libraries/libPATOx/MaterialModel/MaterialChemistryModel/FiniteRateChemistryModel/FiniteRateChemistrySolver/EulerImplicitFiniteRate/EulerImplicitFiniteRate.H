/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::EulerImplicitFiniteRate

Description
    An Euler implicit solver for chemistry

SourceFiles
    EulerImplicitFiniteRate.C

\*---------------------------------------------------------------------------*/

#ifndef EulerImplicitFiniteRate_H
#define EulerImplicitFiniteRate_H

#include "FiniteRateChemistrySolver.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class Type>
class simpleMatrix;

/*---------------------------------------------------------------------------*\
                       Class EulerImplicitFiniteRate Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class EulerImplicitFiniteRate
:
public FiniteRateChemistrySolver<ChemistryModel>
{
  // Private Data

  //- Coefficients dictionary
  dictionary coeffsDict_;


  // Model constants

  //- Chemistry timescale
  scalar cTauChem_;

  //- Equilibrium rate limiter flag (on/off)
  Switch eqRateLimiter_;

  // Solver data
  mutable scalarField cTp_;


 public:

  //- Runtime type information
  TypeName("EulerImplicitFiniteRate");


  // Constructors

  //- Construct from thermo
  EulerImplicitFiniteRate(typename ChemistryModel::reactionThermo& thermo);


  //- Destructor
  virtual ~EulerImplicitFiniteRate();


  // Member Functions

  void updateRRInReactionI
  (
      const label index,
      const scalar pr,
      const scalar pf,
      const scalar corr,
      const label lRef,
      const label rRef,
      const scalar p,
      const scalar T,
      simpleMatrix<scalar>& RR
  ) const;

  //- Update the concentrations and return the chemical time
  virtual void solve
  (
      scalarField& c,
      scalar& T,
      scalar& p,
      scalar& deltaT,
      scalar& subDeltaT
  ) const;

  //- Update the concentrations and return the chemical time, do not update temperature (PATO)
  virtual void solveFixedT
  (
      scalarField &c,
      scalar& T,
      scalar& p,
      scalar& deltaT,
      scalar& subDeltaT
  ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "EulerImplicitFiniteRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
