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
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"MaterialModel",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
materialDict_(
    IOobject
    (
        dictName+"Properties",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
debug_(materialDict_.lookupOrDefault<Switch>("debug", false)),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
materialChemistryModel_(meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh,dictName,"MaterialChemistry")),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,dictName,"MaterialProperties")),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
massModel_(meshLookupOrConstructModel<simpleMassModel>(mesh,dictName,"Mass")),
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh,dictName,"Energy")),
solidMechanicsModel_(meshLookupOrConstructModel<simpleSolidMechanicsModel>(mesh,dictName,"SolidMechanics")),
volumeAblationModel_(meshLookupOrConstructModel<simpleVolumeAblationModel>(mesh,dictName,"VolumeAblation")),
ioModel_(meshLookupOrConstructModel<simpleIOModel>(mesh,dictName,"IO")),
timeControlModel_(meshLookupOrConstructModel<simpleTimeControlModel>(mesh,dictName,"TimeControl"))
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
  }
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

