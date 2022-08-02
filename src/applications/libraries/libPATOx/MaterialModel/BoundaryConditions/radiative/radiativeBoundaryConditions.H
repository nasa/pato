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
    Foam::radiativeBoundaryConditions

Description
    BoundaryConditions using 'radiative'

SourceFiles
    radiativeBoundaryConditions.H
    radiativeBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef radiativeBoundaryConditions_H
#define radiativeBoundaryConditions_H


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
                     Class radiativeBoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

class radiativeBoundaryConditions
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

  //- Convective heat flux at the walls [W/m2]
  volScalarField& qConv_;

  //- Conductive heat flux at the walls [W/m2]
  volScalarField& qCond_; // updated in the surfaceEnergyBalancePatchScalarField

  //- heat transfer coefficient at the walls [kg/m3/s]
  volScalarField& rhoeUeCH;

  //- recovery enthalpy [J/kg]
  volScalarField& h_r;

  //- background temperature [K]
  volScalarField& Tbackground;

  //- scaling factor [-]
  volScalarField& lambda;

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

  //  Stefan–Boltzmann constant [W/m2/K4]
  const dimensionedScalar sigmaSB;

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
  radiativeBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& dictName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~radiativeBoundaryConditions();

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
