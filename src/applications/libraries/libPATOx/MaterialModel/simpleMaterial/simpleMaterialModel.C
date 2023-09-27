/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "simpleMaterialModel.H"

//* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleMaterialModel, 0);
defineRunTimeSelectionTable(simpleMaterialModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMaterialModel::simpleMaterialModel
(
    const fvMesh& mesh,
    const word& regionName
):
IOdictionary
(
    IOobject
    (
        regionName+"MaterialModel",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
#if defined(FOAM_EXTEND)
foam_extend_mesh_(mesh_.time().db().lookupObject<simpleMaterialsModel>("MaterialsModel").get_foam_extend_mesh(regionName)),
#endif
materialDict_(
    IOobject
    (
        regionName+"Properties",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
debug_(materialDict_.lookupOrDefault<Switch>("debug", false)),
regionName_(regionName),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
cellMotionU_ptr(dynamicMesh_?&const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>("cellMotionU")):nullptr),
materialChemistryModel_(meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh,regionName,simpleMaterialChemistryModel::modelName)),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,regionName,simplePyrolysisModel::modelName)),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,regionName,simpleMaterialPropertiesModel::modelName)),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,regionName,simpleGasPropertiesModel::modelName)),
massModel_(meshLookupOrConstructModel<simpleMassModel>(mesh,regionName,simpleMassModel::modelName)),
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh,regionName,simpleEnergyModel::modelName)),
solidMechanicsModel_(meshLookupOrConstructModel<simpleSolidMechanicsModel>(mesh,regionName,simpleSolidMechanicsModel::modelName)),
volumeAblationModel_(meshLookupOrConstructModel<simpleVolumeAblationModel>(mesh,regionName,simpleVolumeAblationModel::modelName)),
ioModel_(meshLookupOrConstructModel<simpleIOModel>(mesh,regionName,simpleIOModel::modelName)),
timeControlModel_(meshLookupOrConstructModel<simpleTimeControlModel>(mesh,regionName,simpleTimeControlModel::modelName))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMaterialModel::~simpleMaterialModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

void Foam::simpleMaterialModel::update()
{
  if(this->dynamicMesh_) {
    if(debug_) {
      Info << "update mesh " << endl;
    }
    dynamicFvMesh& mesh = (dynamicFvMesh&)(this->mesh_);
    mesh.update();
    const Time& runTime = mesh.time();
#include "volContinuity.H" // checks the topology of the mesh after the motion
    volVectorField& cellMotionU_ = *cellMotionU_ptr;
    forAll(cellMotionU_, cellI) {
      cellMotionU_[cellI] = vector::zero;
    }
    forAll(cellMotionU_.boundaryField(), patchI) {
      forAll(cellMotionU_.boundaryField()[patchI], faceI) {
        cellMotionU_.boundaryFieldRef()[patchI][faceI] = vector::zero;
      }
    }
  }
#if defined(FOAM_EXTEND)
  Foam_extend_::scalar time_value = mesh_.time().value(); // Foam
  const_cast<Foam_extend_::Time&>(foam_extend_mesh_.time()).setTime(time_value,0); // Foam_extend
  if (this->dynamicMesh_) {
    forAll(foam_extend_mesh_.points() , i) {
      forAll(foam_extend_mesh_.points()[i], j) {
        double& extend_point_j = const_cast<double&>(foam_extend_mesh_.points()[i][j]); // Foam_extend
        extend_point_j=mesh_.points()[i][j];
      }
    }
  }
#endif
  if(debug_) {
    Info << "update materialPropertiesModel " << endl;
  }
  materialPropertiesModel_.update();
  if(debug_) {
    Info << "update gasPropertiesModel " << endl;
  }
  gasPropertiesModel_.update();
  if(debug_) {
    Info << "update MaterialChemistryModel_ " << endl;
  }
  materialChemistryModel_.update();
  if(debug_) {
    Info << "update pyrolysisModel_ " << endl;
  }
  pyrolysisModel_.update();
  if(debug_) {
    Info << "update heterogeneousModel_ " << endl;
  }
  volumeAblationModel_.update();
  if(debug_) {
    Info << "update massModel_ " << endl;
  }
  massModel_.update();
  if(debug_) {
    Info << "update energyModel_ " << endl;
  }
  energyModel_.update();
  if(debug_) {
    Info << "update solidMechanicsModel_ " << endl;
  }
  solidMechanicsModel_.update();
  if(debug_) {
    Info << "update ioModel_ " << endl;
  }
  ioModel_.update();
  if(debug_) {
    Info << "update timeControlModel_ " << endl;
  }
  timeControlModel_.update();
  if(debug_) {
    Info << "end --- Foam::simpleMaterialModel::update() " << endl;
  }
  updateMinDeltaTw_ = timeControlModel_.updateMinDeltaTw(); // NAmrofel 26/11/2019
}

// ************************************************************************* //

