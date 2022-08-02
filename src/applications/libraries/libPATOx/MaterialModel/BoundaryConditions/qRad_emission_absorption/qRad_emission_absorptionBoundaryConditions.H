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
    along with OpenFOAM.  If radiativet, see <http://www.gnu.org/licenses/>.

Class
    Foam::qRad_emission_absorptionBoundaryConditions

Description
    BoundaryConditions using 'radiative'

SourceFiles
    qRad_emission_absorptionBoundaryConditions.H
    qRad_emission_absorptionBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef qRad_emission_absorptionBoundaryConditions_H
#define qRad_emission_absorptionBoundaryConditions_H


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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class qRad_emission_absorptionBoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

class qRad_emission_absorptionBoundaryConditions
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

  //- Re-radiative heat flux at the walls [W/m2]
  volScalarField& qRadEmission_;

  //- Radiative heat flux from the environment at the walls [W/m2]
  volScalarField& qRadAbsorption_;

  //- Conductive heat flux at the walls [W/m2]
  volScalarField& qCond_; // updated in the surfaceEnergyBalancePatchScalarField

  //- background temperature [K]
  volScalarField& Tbackground;

  //- radiative heat flux from the environment at the walls [W/m3] (qRad)
  volScalarField& qRad;

  //- Solid emissivity [-]
  volScalarField& emissivity_;

  //- Solid absorptivity [-]
  volScalarField& absorptivity_;

  //- Solid conductivity [W/m/K]
  const volTensorField& k_;

  //- Temperature (K)
  volScalarField& T_;

  //- Pressure (Pa)
  volScalarField& p_;

  //  Stefan–Boltzmann constant [W/m2/K4]
  const dimensionedScalar sigmaSB;

  //- factor for radiative heating [-]
  scalar qRad_factor_;

  //- Debug switch
  Switch debug_;

  //- Needed fields in boundaryMapping
  wordList neededFields_;

  //- BoundaryMapping object
  autoPtr<simpleBoundaryMappingModel> boundaryMapping_;
  simpleBoundaryMappingModel * boundaryMapping_ptr;

 public:

  // Constructors

  //- Construct
  qRad_emission_absorptionBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& dictName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~qRad_emission_absorptionBoundaryConditions();

  // Member Functions
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
