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
    along with OpenFOAM.  If HeatFluxt, see <http://www.gnu.org/licenses/>.

Class
    Foam::HeatFluxBoundaryConditions

Description
    BoundaryConditions using 'HeatFlux'

SourceFiles
    HeatFluxBoundaryConditions.H
    HeatFluxBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef HeatFluxBoundaryConditions_H
#define HeatFluxBoundaryConditions_H


#undef Log // conflict between OpenFoam and Mutation++ on the alias "Log" -> undefined here from OpenFoam.
#include <mutation++/mutation++.h>
#include "fvCFD.H"
#include "mathFunctions.H"
#include "IOmanip.H"
#include "IOFunctions.H"
#include "simpleBoundaryMappingModel.H"
#include "PATOx.H"
#include "fluxFactor.H"
#include "boundaryMappingFvPatchScalarField.H"
#include "optionsBC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class HeatFluxBoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

class HeatFluxBoundaryConditions
{
 protected:

  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary name
  const word regionName_;

  //- Energy model
  simpleEnergyModel& energyModel_;

  //- Initialized flag
  Switch initialized_;

  //- volScalarFields
  PMap<volScalarField> scalarFields_;

  //- volVectorFields
  PMap<volVectorField> vectorFields_;

  //- volTensorFields
  PMap<volTensorField> tensorFields_;

  //- Current patch ID
  const label currentPatchID_;

  //- Dictionary name
  const dictionary dict_;

  //- True if mesh is a dynamic mesh
  Switch dynamicMesh_;

  //  Stefan–Boltzmann constant [W/m2/K4]
  const dimensionedScalar sigmaSB;

  // Component of dynamic mesh
  word component_;

  //- Debug switch
  Switch debug_;

  //- Needed fields in boundaryMapping
  wordList neededFields_;

  //- BoundaryMapping object
  autoPtr<simpleBoundaryMappingModel> boundaryMapping_;

  simpleBoundaryMappingModel * boundaryMapping_ptr;

  //- Options BC
  optionsBC optionsBC_;

 public:

  // Constructors

  //- Construct
  HeatFluxBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& regionName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~HeatFluxBoundaryConditions();

  // Member Functions

  //- Initialize
  void initialize();

  //- Update
  virtual void update();
  inline void updateTemperatureBC();
  inline void updateBoundaryMapping();

  //- Fields
  inline wordList neededFields();

  //- write in os
  void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
