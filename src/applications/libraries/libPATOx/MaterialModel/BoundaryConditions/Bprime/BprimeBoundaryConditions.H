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
    Foam::BprimeBoundaryConditions

Description
    BoundaryConditions using 'Bprime'

SourceFiles
    BprimeBoundaryConditions.H
    BprimeBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef BprimeBoundaryConditions_H
#define BprimeBoundaryConditions_H


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
#include "optionsBC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BprimeBoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

class BprimeBoundaryConditions
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

  // Bprime mixture from Mutation++

  //- Name of the mixture
  word mixtureMutationBprime;

  //- Environment directory
  fileName environmentDirectory;

  //- Environment dictionary
  IOdictionary environmentDictionary;

  //- Inputs for the mole fractions of the gas composition inside the material. Not necessary if there are Zx[...] in constantProperties
  List<Tuple2<word,scalar> > moleFractionGasInMaterial_;

  //- Switch to write Bprime table
  Switch writeBprimeTable_;

  //- Minimum BprimeC
  scalar minBc;

  //- Maximum BprimeC
  scalar maxBc;

  //- Switch to read Bprime table
  Switch readBprimeTable_;

  //- Object BprimeTable
  BprimeTable * BprimeT_ptr;

  //- Pyrolysis gas mass flow rate at the walls [kg/m2/s]
  volScalarField& mDotGw_;

  //- Char ablation rate at the walls [kg/m2/s]
  volScalarField& mDotCw_;

  //- Dimensionless char ablation rate at the walls
  volScalarField& BprimeCw_;

  //- Dimensionless pyrolysis gas blowing rate at the walls
  volScalarField& BprimeGw_;

  //- Re-radiative heat flux at the walls [W/m2]
  volScalarField& qRadEmission_;

  //- Radiative heat flux from the environment at the walls [W/m2]
  volScalarField& qRadAbsorption_;

  //- Advective heat flux at the walls due to the pyrolysis gas [W/m2]
  volScalarField& qAdvPyro_;

  //- Advective heat flux at the walls due to the char [W/m2]
  volScalarField& qAdvChar_;

  //- Convective heat flux at the walls [W/m2]
  volScalarField& qConv_;

  //- Conductive heat flux at the walls [W/m2]
  volScalarField& qCond_; // updated in the surfaceEnergyBalancePatchScalarField

  //- Recession [m]
  volScalarField& recession_;

  //- heat transfer coefficient at the walls [kg/m2/s]
  volScalarField& rhoeUeCH;

  //- blowing correction due to the blockage induced by the pyrolysis and ablation gas blowing [-]
  volScalarField& blowingCorrection;

  //- wall enthalpy [J/kg]
  volScalarField& h_w;

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

  //- Pyrolysis gas mass flow rate [kg/m2/s]
  const volVectorField& mDotG_;

  //- Pyrolysis gas mass flux flow rate interpolated on the faces [kg/m2/s]
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

  //- Solid density [-]
  const volScalarField& rho_s_;

  //- Pressure [Pa]
  volScalarField& p_;

  //- Temperature [K]
  volScalarField& T_;

  //- Material dict
  const IOdictionary materialDict_;

  //- material properties directory
  const fileName materialPropertiesDirectory;

  //- constantProperties dictionary
  IOdictionary constantPropertiesDictionary;

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

  //- Mass fractions of the environment composition [-]
  double * Yke_ref;

  //- Species mass fraction of the gas composition inside the material [-]
  double * Ykg_ref;

  //- Options BC
  optionsBC optionsBC_;

 public:
  //- Mutation++ Mixture
  Mutation::Mixture * mixBprime_ptr;

  //- Number of elements in the mixture
  int ne_Bprime;

  //- Number of species in the mixture
  int ns_Bprime;

  //- volScalarField names of the species mass fraction of the gas composition inside the material [-]
  wordList Ykg_names;

  //- volScalarField names of the species mole fraction of the gas composition at the wall [-]
  wordList Xw_names;

  //- If recession flag = false, there is no recession
  Switch recessionFlag;

  //- Bprime table filename to read. If not there, compute B' on-the-fly
  fileName nameBprimeFile_;

  // Constructors

  //- Construct
  BprimeBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& dictName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~BprimeBoundaryConditions();

  // Member Functions
  //- Bprime
  inline void initBprimeMixture();
  inline void writeBprimeTable();
  inline void readBprimeTable();
  virtual void computeSurfaceMassBalance(scalarList Yke, scalarList Ykg, scalar pw, scalar Tw, scalar Bg, scalar& Bc_, scalar& hw_, scalarList& p_Xw_);

  //- Environment
  inline void readEnvironmentTable();

  //- Update
  virtual void update();
  inline void updateBprimeBC();
  inline void updateTemperatureBC();
  inline void updateMotionBC();
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
