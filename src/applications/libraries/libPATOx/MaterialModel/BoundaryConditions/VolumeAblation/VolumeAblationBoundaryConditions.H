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
    Foam::VolumeAblationBoundaryConditions

Description
    BoundaryConditions using 'Bprime'

SourceFiles
    VolumeAblationBoundaryConditions.H
    VolumeAblationBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef VolumeAblationBoundaryConditions_H
#define VolumeAblationBoundaryConditions_H


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
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "simpleGasPropertiesModel.H"
#include "simpleMaterialChemistryModel.H"
#include "simpleVolumeAblationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class VolumeAblationBoundaryConditions Declarations
\*---------------------------------------------------------------------------*/

class VolumeAblationBoundaryConditions
{
 protected:

  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary name
  const word regionName_;

  //- Current patch ID
  const label currentPatchID_;

  //- Dictionary name
  const dictionary dict_;

  //- Boundary name
  const word boundaryName_;

  //- Energy Model
  simpleEnergyModel& energyModel_;

  //- Initialized flag
  Switch initialized_;

  //- volScalarFields
  PMap<volScalarField> scalarFields_;

  //- volVectorFields
  PMap<volVectorField> vectorFields_;

  //- volTensorFields
  PMap<volTensorField> tensorFields_;

  //- True if mesh is a dynamic mesh
  Switch dynamicMesh_;

  //- pressure and temperature pointer for Mutation++
  double * p_pT;

  //- Physics based erosion model
  const int physicsBasedErosionModel_;

  //- Environment directory
  fileName environmentDirectory;

  //- Environment dictionary
  IOdictionary environmentDictionary_;

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

  //  Stefan–Boltzmann constant [W/m2/K4]
  const dimensionedScalar sigmaSB;

  // Component of dynamic mesh
  word component_;

  // cell motion velocity of the dynamic mesh
  volVectorField * cellMotionU_ptr;

  // cell motion velocity of the dynamic mesh
  volScalarField * cellMotionUxyz_ptr;

  //- initial position of the moving surface [m]
  volVectorField initialPosition_;

  //- Debug switch
  Switch debug_;

  //- Needed fields in boundaryMapping
  wordList neededFields_;

  //- BoundaryMapping object
  autoPtr<simpleBoundaryMappingModel> boundaryMapping_;
  simpleBoundaryMappingModel * boundaryMapping_ptr;

  //- Switch that avoids deleting empty pointers if the destructor is called twice
  static Switch ptr_deleted;

 public:

  //- Mass fractions of the environment (boundary layer edge) specie composition [-]
  double * Yie_ref;

  // Constructors

  //- Construct
  VolumeAblationBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& regionName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~VolumeAblationBoundaryConditions();

  // Member Functions

  //- Boundary name
  inline word boundaryName();

  //- Initialize
  inline void initialize();

  //- Update
  virtual void update();
  inline void updateTemperatureBC();
  inline void updateMotionBC();
  inline void updateBoundaryMapping();

  //- Fields
  inline wordList neededFields();

  //- write in os
  void write(Ostream& os) const;

  // return members

  //- heat transfer coefficient at the walls [kg/m3/s]
  inline volScalarField& rhoeUeCH() ;

  //- blowing correction due to the blockage induced by the pyrolysis and ablation gas blowing [-]
  inline volScalarField& blowingCorrection();

  //- Environment dictionary
  inline IOdictionary& environmentDictionary();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VolumeAblationBoundaryConditionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
