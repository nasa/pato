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
    along with OpenFOAM.  If PurePyrolysist, see <http://www.gnu.org/licenses/>.

Class
    Foam::PurePyrolysisIOModel

Description

SourceFiles
    PurePyrolysisIOModel.H
    PurePyrolysisIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef PurePyrolysisIOModel_H
#define PurePyrolysisIOModel_H

#include "simpleIOModel.H"
#include "simpleMaterialChemistryModel.H"
#include "simplePyrolysisModel.H"
#include "simpleMaterialPropertiesModel.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PurePyrolysisIOModel Declaration
\*---------------------------------------------------------------------------*/

class PurePyrolysisIOModel
:
public simpleIOModel
{

 protected:

  //- List of wordList null
  const List<wordList> nullListList_;

  //- runTimeModifiable in controlDict
  Switch runTimeModifiable_;

  //- writeInterval in controlDict
  scalar writeInterval_;

  //- writeControl in controlDict
  word writeControl_;

  //- Temperature [K]
  volScalarField& Ta_;

  //- MaterialChemistry Model
  simpleMaterialChemistryModel& MaterialChemistryModel_;

  //- Elements names
  const wordList& elementNames_;

  //- Species names
  const wordList& speciesNames_;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Pyrolysis reaction rates per specie/element [kg/m^3/s]
  const PtrList<volScalarField>& pi_;

  //- Total pyrolysis reaction rates [kg/m^3/s]
  const volScalarField& piTotal_;

  //- Material Properties Model
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- Solid mass density [kg/m3]
  const volScalarField& rho_s_;

  //- Initial solid mass density [kg/m3]
  const dimensionedScalar initial_rho_s_;

  // Output Model: READ
  // List of the fields to write in the time folders
  const wordList& outputList_;
  //- List of the files to read
  const List<fileName>& readFilesList_;
  //- List of the data in the read files
  const List<List<scalarList> >& filesData_;
  //- Data from the readFiles: Species pyrolysis gas production rate [kg/m3/s]
  PtrList<List<scalarList> > speciesData_;
  //- Data from the readFiles: Ratio of the solid density (rho_s/rho_sv) [-]
  List<scalarList> rhoRatioData_;
  //- Data from the readFiles: Temperature [K]
  List<scalarList> temperatureData_;
  //- Data from the readFiles: Normalized total pyrolysis gas production rate (piTot/rho_sv) [1/s]
  List<scalarList> normalizedPiTotalData_;
  //- File names of the "readFiles" for the rhoRatio, T and normPiTot
  const word rhoRatioFileName_;
  const word temperatureFileName_;
  const word normalizedPiTotalFileName_;

  // Output Model: COMPUTE
  //- Ratio of the solid density (rho_s/rho_sv) [-]
  volScalarField rhoRatio_;
  //-  Normalized total pyrolysis gas production rate (piTot/rho_sv) [1/s]
  volScalarField normalizedPiTotal_;
  //- Species mass fractions y_i [-]
  PtrList<volScalarField>& massFractions_;
  //- Normalized species pyrolysis gas production rates (y_i*piTot/rho_sv) [1/s]
  PtrList<volScalarField> normalizedPi_;
  //- Switch to initialize the normalizePi_
  Switch initNormalizedPi_;

  // OutputModel: COMPARE
  // list of the (file fields) to compare
  const List<wordList> compareList_;
  // list of the fields to compare
  wordList compareFields_;
  // list of the data file to compare
  wordList compareData_;
  // First time in the compare files
  scalar compareFirstTime_;
  // Dakota results file name
  const fileName dakotaResultsFileName_;
  // Dakota responses (objective or calibration)
  const word dakotaResponses_;
  // Dakota results
  scalarList dakotaResults_;
  // index of the species (comparison species order from compareList and PATO)
  List<int> indexSpeciesCompareData_;
  // number of pre-run
  int preRunNumber_;

  // OutputModel: INITIALIZE
  //- Switch to initialize the basic output model
  const bool initOutput_;

 public:

  //- Runtime type information
  TypeName("PurePyrolysis");

  // Constructors

  //- Construct from mesh and phase name
  PurePyrolysisIOModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PurePyrolysisIOModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Write output
  void writeOutput();

  //- Read the input files
  void readInput();

  //- Initialize normalized pi field
  Switch initNormalizedPi();

  //- Update the sampling
  void updateSampling();

  //- Update the delta time
  void updateDeltaTime(const scalarList& times);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
