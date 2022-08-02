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
  const word dictName_;

  //- Current patch ID
  const label currentPatchID_;

  //- Dictionary name
  const dictionary dict_;

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

  //- Re-radiative heat flux at the walls [W/m2]
  volScalarField& qRadEmission_;

  //- Radiative heat flux from the environment at the walls [W/m2]
  volScalarField& qRadAbsorption_;

  //- Convective heat flux at the walls [W/m2]
  volScalarField& qConv_;

  //- Conductive heat flux at the walls [W/m2]
  volScalarField& qCond_; // updated in the surfaceEnergyBalancePatchScalarField

  //- Recession [m]
  volScalarField& recession_;

  //- heat transfer coefficient at the walls [kg/m3/s]
  volScalarField& rhoeUeCH_;

  //- blowing correction due to the blockage induced by the pyrolysis and ablation gas blowing [-]
  volScalarField& blowingCorrection_;

  //- wall enthalpy using the gaseous species composition in the porous material [J/kg]
  volScalarField& h_w;

  //- wall enthalpy using the gaseous species composition of the boundayr layer edge [J/kg]
  volScalarField& h_ew;

  //- recovery enthalpy [J/kg]
  volScalarField& h_r;

  //- background temperature [K]
  volScalarField& Tbackground;

  //- scaling factor [-]
  volScalarField& lambda;

  //- radiative heat flux from the environment at the walls [W/m3] (qRad)
  volScalarField& qRad;

  //- Heat On flag
  volScalarField& heatOn;

  //- Pyrolysis gas mass flow rate [kg/m3/s]
  const volVectorField& mDotG_;

  //- Pyrolysis gas mass flux flow rate interpolated on the faces [kg/m3/s]
  const surfaceVectorField& mDotGFace_;

  //- Pyrolysis gas enthalpy [J/kg]
  const volScalarField& h_g_;

  //- Char enthalpy [J/kg]
  const volScalarField& h_c_;

  //- Solid emissivity [-]
  const volScalarField& emissivity_;

  //- Solid absorptivity [-]
  const volScalarField& absorptivity_;

  //- Solid conductivity [W/m/K]
  const volTensorField& k_;

  //- Solid fiber radius [m]
  const volScalarField& rT;

  //- Heterogeneous rate [kg/m3/s]
  const volScalarField& omegaHeterogeneousRate_;

  //- Solid density [kg/m3]
  const volScalarField& rho_s_;

  //- Pressure (Pa)
  volScalarField& p_;

  //- Temperature (K)
  volScalarField& T_;

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

 public:
  //- Mass fractions of the environment (boundary layer edge) specie composition [-]
  double * Yie_ref;

  // Constructors

  //- Construct
  VolumeAblationBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& dictName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~VolumeAblationBoundaryConditions();

  // Member Functions
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
  inline volScalarField& rhoeUeCH() const;

  //- blowing correction due to the blockage induced by the pyrolysis and ablation gas blowing [-]
  inline volScalarField& blowingCorrection() const;

  inline IOdictionary& environmentDictionary() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VolumeAblationBoundaryConditionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
