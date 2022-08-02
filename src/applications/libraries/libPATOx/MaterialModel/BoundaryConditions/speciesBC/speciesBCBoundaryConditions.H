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
    along with OpenFOAM.  If Bprimet, see <http://www.gnu.org/licenses/>.

Class
    Foam::speciesBCBoundaryConditions

Description
    BoundaryConditions using 'Bprime'

SourceFiles
    speciesBCBoundaryConditions.H
    speciesBCBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef speciesBCBoundaryConditions_H
#define speciesBCBoundaryConditions_H


#undef Log // conflict between OpenFoam and Mutation++ on the alias "Log" -> undefined here from OpenFoam.
#include <mutation++/mutation++.h>
#include "fvCFD.H"
#include "BprimeTable.H"
#include "mathFunctions.H"
#include "IOmanip.H"
#include "IOFunctions.H"
#include "simpleBoundaryMappingModel.H"
#include "PATOx.H"
#include "fluxFactor.H"
#include "boundaryMappingFvPatchScalarField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class speciesBCBoundaryConditions Declarations
\*---------------------------------------------------------------------------*/

class speciesBCBoundaryConditions
{
 protected:

  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary name
  const word dictName_;

  //- Current patch ID
  const label currentPatchID_;

  //- Dictionary name
  const dictionary dict_;

  //- True if mesh is a dynamic mesh
  Switch dynamicMesh_;

  //- Environment directory
  fileName environmentDirectory;

  //- Environment dictionary
  IOdictionary environmentDictionary_;

  //- heat transfer coefficient at the walls [kg/m3/s]
  volScalarField& rhoeUeCH_;

  //- blowing correction due to the blockage induced by the pyrolysis and ablation gas blowing [-]
  volScalarField& blowingCorrection_;

  //- Material dict
  const IOdictionary materialDict_;

  //- Mutation++ Mixture name
  word mixtureMutation;

  //- Mutation++ Mixture Options
  Mutation::MixtureOptions * optEq_;

  //- Mutation++ Mixture
  Mutation::Mixture * mix_;

  //- Numer of species in the Mutation++ mixture
  int ns_mix;

  //- material properties directory
  const fileName materialPropertiesDirectory;

  //- constantProperties dictionary
  IOdictionary constantPropertiesDictionary;

  //- Failure fiber radius [m]
  dimensionedScalar rfFail;

  //- Debug switch
  Switch debug_;

 public:
  //- Mass fractions of the environment (boundary layer edge) specie composition [-]
  double * Yie_ref;

  // Constructors

  //- Construct
  speciesBCBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& dictName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~speciesBCBoundaryConditions();

  // Member Functions
  //- Update
  virtual void update();

  //- write in os
  void write(Ostream& os) const;

  // return members

  //- return the environment dictionary
  inline IOdictionary& environmentDictionary() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "speciesBCBoundaryConditionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
