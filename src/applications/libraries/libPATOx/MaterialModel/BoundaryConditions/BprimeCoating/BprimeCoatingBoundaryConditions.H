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
    Foam::BprimeCoatingBoundaryConditions

Description
    BoundaryConditions using 'Bprime'

SourceFiles
    BprimeCoatingBoundaryConditions.H
    BprimeCoatingBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef BprimeCoatingBoundaryConditions_H
#define BprimeCoatingBoundaryConditions_H


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
                     Class BprimeCoatingBoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

class BprimeCoatingBoundaryConditions
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

  //- Material dictionary
  const IOdictionary& materialDict_;

  //- True if mesh is a dynamic mesh
  Switch dynamicMesh_;

  //- Bprime table filename to read. If not there, compute B' on-the-fly
  fileName nameBprimeFile_;

  //- Object BprimeTable
  BprimeTable * BprimeT_ptr;

  //- Object coating BprimeTable
  BprimeTable * coatingBprimeT_ptr;

  //  Stefan–Boltzmann constant [W/m2/K4]
  const dimensionedScalar sigmaSB;

  // Component of dynamic mesh
  word component_;

  // cell motion velocity of the dynamic mesh
  volVectorField * cellMotionU_ptr;

  // cell motion velocity of the dynamic mesh
  volScalarField * cellMotionUxyz_ptr;

  //- initial position of the moving surface [m]
  volVectorField& initialPosition_;

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

  //- Coating recession [m]
  volScalarField& coatingRecession_;

  //- Coating height
  scalar coatingHeight_;

  //- Coating table file name
  fileName coatingTableFileName_;

  //- Coating table
  List<scalarList> coatingTable_;

  //- Coating Bprime file name
  fileName coatingBprimeFileName_;

  //- MaterialProperties type
  word typeNameMaterialPropertiesModel;

  //- Option BC for B'c and h_w
  optionsBC option_BprimeC_hw_;

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

  // Constructors

  //- Construct
  BprimeCoatingBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& regionName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~BprimeCoatingBoundaryConditions();

  // Member Functions
  //- Bprime
  inline void initBprimeMixture();
  inline RectangularMatrix<scalar> readBprimeTable(fileName BprimeFile);
  virtual void computeSurfaceMassBalance(label faceI, scalarList Yke, scalarList Ykg, scalar pw, scalar Tw, scalar Bg, scalar& Bc_, scalar& hw_, scalarList& p_Xw_);

  //- Environment
  inline void readEnvironmentTable();

  //- Initialize
  void initialize();

  //- Update
  virtual void update();
  inline void updateCoatingBC();
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
