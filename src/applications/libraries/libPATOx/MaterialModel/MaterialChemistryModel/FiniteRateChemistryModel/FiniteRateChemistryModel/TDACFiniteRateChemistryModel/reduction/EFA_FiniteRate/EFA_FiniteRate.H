/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::FiniteRateChemistryReductionMethods::EFA_FiniteRate

Description

SourceFiles
    EFA_FiniteRate.C

\*---------------------------------------------------------------------------*/

#ifndef EFA_FiniteRate_H
#define EFA_FiniteRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FiniteRateChemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                           Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class EFA_FiniteRate
:
public FiniteRateChemistryReductionMethod<CompType, ThermoType>
{
  // Private Data

  List<label> sC_,sH_,sO_,sN_;
  scalar  sortPart_;

 public:

  //- Runtime type information
  TypeName("EFA_FiniteRate");


  // Constructors

  //- Construct from components
  EFA_FiniteRate
  (
      const IOdictionary& dict,
      TDACFiniteRateChemistryModel<CompType, ThermoType>& chemistry
  );


  //- Destructor
  virtual ~EFA_FiniteRate();


  // Member Functions

  //- Reduce the mechanism
  virtual void reduceMechanism
  (
      const scalarField &c,
      const scalar T,
      const scalar p
  );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FiniteRateChemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "EFA_FiniteRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
