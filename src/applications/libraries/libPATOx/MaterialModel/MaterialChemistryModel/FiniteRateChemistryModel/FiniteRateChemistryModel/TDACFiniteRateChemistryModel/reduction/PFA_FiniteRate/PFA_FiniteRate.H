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
    Foam::FiniteRateChemistryReductionMethods::PFA_FiniteRate

Description
    Path flux analysis

SourceFiles
    PFA_FiniteRate.C

\*---------------------------------------------------------------------------*/

#ifndef PFA_FiniteRate_H
#define PFA_FiniteRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FiniteRateChemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                           Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class PFA_FiniteRate
:
public FiniteRateChemistryReductionMethod<CompType, ThermoType>
{
  // Private Data

  //- List of label for the search initiating set
  labelList searchInitSet_;


 public:

  //- Runtime type information
  TypeName("PFA_FiniteRate");


  // Constructors

  //- Construct from components
  PFA_FiniteRate
  (
      const IOdictionary& dict,
      TDACFiniteRateChemistryModel<CompType, ThermoType>& chemistry
  );


  //- Destructor
  virtual ~PFA_FiniteRate();


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
#include "PFA_FiniteRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
