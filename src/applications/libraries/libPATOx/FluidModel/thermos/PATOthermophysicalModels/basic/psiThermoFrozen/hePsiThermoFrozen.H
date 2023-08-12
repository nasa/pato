/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::hePsiThermoFrozen

Description
    Energy for a mixture based on compressibility

SourceFiles
    hePsiThermoFrozen.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiThermoFrozen_H
#define hePsiThermoFrozen_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hePsiThermoFrozen Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class hePsiThermoFrozen
:
public heThermo<BasicPsiThermo, MixtureType>
{
  // Private Member Functions

  //- Calculate the thermo variables
  void calculate();

  //- Construct as copy (not implemented)
  hePsiThermoFrozen(const hePsiThermoFrozen<BasicPsiThermo, MixtureType>&);

 public:

  //- Runtime type information
  TypeName("hePsiThermoFrozen");


  // Constructors

  //- Construct from mesh and phase name
  hePsiThermoFrozen
  (
      const fvMesh&,
      const word& phaseName
  );


  //- Destructor
  virtual ~hePsiThermoFrozen();


  // Member functions

  //- Update properties
  virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "hePsiThermoFrozen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
