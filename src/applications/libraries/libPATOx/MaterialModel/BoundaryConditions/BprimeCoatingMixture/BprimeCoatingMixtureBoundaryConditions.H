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
    along with OpenFOAM.  If BprimeCoatingMixturet, see <http://www.gnu.org/licenses/>.

Class
    Foam::BprimeCoatingMixtureBoundaryConditions

Description
    BoundaryConditions using 'BprimeCoatingMixture'

SourceFiles
    BprimeCoatingMixtureBoundaryConditions.H
    BprimeCoatingMixtureBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef BprimeCoatingMixtureBoundaryConditions_H
#define BprimeCoatingMixtureBoundaryConditions_H


#undef Log // conflict between OpenFoam and Mutation++ on the alias "Log" -> undefined here from OpenFoam.
#include <mutation++/mutation++.h>
#include "fvCFD.H"
#include "BprimeCoatingMixtureTable.H"
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
                     Class BprimeCoatingMixtureBoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

class BprimeCoatingMixtureBoundaryConditions
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

  //- Material dict
  const IOdictionary& materialDict_;

  //- True if mesh is a dynamic mesh
  Switch dynamicMesh_;

  //- Name of the mixture
  word mixtureName;

  //- Environment directory
  fileName environmentDirectory;

  //- Environment dictionary
  IOdictionary environmentDictionary;

  //- Minimum BprimeC
  scalar minBc;

  //- Maximum BprimeC
  scalar maxBc;

  //  Stefan–Boltzmann constant [W/m2/K4]
  const dimensionedScalar sigmaSB;

  // cell motion velocity of the dynamic mesh
  volVectorField * cellMotionU_ptr;

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

  //- Element composition on the surface
  List<Tuple2<string,double>> elemOnSurface;

  //- Mode of the surface compostion
  int mode_surf;

 public:
  //- Mutation++ Mixture
  Mutation::Mixture * mixture_ptr;

  //- Number of elements in the mixture
  int ne_mix;

  //- Number of species in the mixture
  int ns_mix;

  //- volScalarField names of the species mass fraction of the gas composition inside the material [-]
  wordList Ykg_names;

  //- volScalarField names of the species mole fraction of the gas composition at the wall [-]
  wordList Xw_names;

  //- If recession flag = false, there is no recession
  Switch recessionFlag;

  // Constructors

  //- Construct
  BprimeCoatingMixtureBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& regionName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~BprimeCoatingMixtureBoundaryConditions();

  // Member Functions
  //- Mixture
  inline void initMixture();
  List<Tuple2<string,double>> computeSiCSiO2SurfComp(const Mutation::Mixture& mix, const double& zC);
  inline Eigen::VectorXd equilibrium(const Mutation::Mixture& mix, const double& p, const double& T, const Eigen::VectorXd& Z_wsg, Eigen::VectorXd& X_sg);
  virtual void computeSurfaceMassBalance(const Mutation::Mixture& mix, const List<Tuple2<string,double>>& elemOnSurface, const Eigen::VectorXd& z_e, const Eigen::VectorXd& z_pg, const double T,
                                         const double p, const double Bg, double &Bc, double &hw, Eigen::VectorXd& X_sg);

  //- Environment
  inline void readEnvironmentTable();

  //- Initialize
  void initialize();

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
