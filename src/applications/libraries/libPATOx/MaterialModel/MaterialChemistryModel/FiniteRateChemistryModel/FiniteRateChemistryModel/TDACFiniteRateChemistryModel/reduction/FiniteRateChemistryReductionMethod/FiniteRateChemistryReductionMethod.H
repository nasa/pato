/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::FiniteRateChemistryReductionMethod

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    FiniteRateChemistryReductionMethod.C

\*---------------------------------------------------------------------------*/

#ifndef FiniteRateChemistryReductionMethod_H
#define FiniteRateChemistryReductionMethod_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CompType, class ThermoType>
class TDACFiniteRateChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class FiniteRateChemistryReductionMethod
{

 protected:

  const IOdictionary& dict_;

  //- Dictionary that store the algorithm data
  const dictionary coeffsDict_;

  //- Is mechanism reduction active?
  Switch active_;

  //- Switch to select performance logging
  Switch log_;

  TDACFiniteRateChemistryModel<CompType, ThermoType>& chemistry_;

  //- List of active species (active = true)
  List<bool> activeSpecies_;

  //- Number of active species
  label NsSimp_;

  //- Number of species
  const label nSpecie_;

  //- Tolerance for the mechanism reduction algorithm
  scalar tolerance_;


 public:

  //- Runtime type information
  TypeName("FiniteRateChemistryReductionMethod");


  // Declare runtime constructor selection table
  declareRunTimeSelectionTable
  (
      autoPtr,
      FiniteRateChemistryReductionMethod,
      dictionary,
      (
          const IOdictionary& dict,
          TDACFiniteRateChemistryModel<CompType, ThermoType>& chemistry
      ),
      (dict, chemistry)
  );


  // Constructors

  //- Construct from components
  FiniteRateChemistryReductionMethod
  (
      const IOdictionary& dict,
      TDACFiniteRateChemistryModel<CompType, ThermoType>& chemistry
  );


  // Selector

  static autoPtr<FiniteRateChemistryReductionMethod<CompType, ThermoType>> New
      (
          const IOdictionary& dict,
          TDACFiniteRateChemistryModel<CompType, ThermoType>& chemistry
      );


  //- Destructor
  virtual ~FiniteRateChemistryReductionMethod();


  // Member Functions

  //- Is mechanism reduction active?
  inline bool active() const;

  //- Is performance data logging enabled?
  inline bool log() const;

  //- Return the active species
  inline const List<bool>& activeSpecies() const;

  //- Return the number of active species
  inline label NsSimp();

  //- Return the initial number of species
  inline label nSpecie();

  //- Return the tolerance
  inline scalar tolerance() const;

  //- Reduce the mechanism
  virtual void reduceMechanism
  (
      const scalarField &c,
      const scalar T,
      const scalar p
  ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FiniteRateChemistryReductionMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "FiniteRateChemistryReductionMethod.C"
#include "FiniteRateChemistryReductionMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
