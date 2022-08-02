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

#include "simpleMaterialsModel.H"

//* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleMaterialsModel, 0);
defineRunTimeSelectionTable(simpleMaterialsModel, Time);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMaterialsModel::simpleMaterialsModel
(
    const Time& runTime
):
IOdictionary
(
    IOobject
    (
        "MaterialsModel",
        runTime.constant(),
        runTime.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
runTime_(runTime),
rp(runTime_),
solidRegionNames_(isFile(changeEnviVar("$FOAM_CASE/constant/regionProperties"))?rp["solid"]:wordList(0))

{
  simpleMinDeltaTw_ = GREAT; // NAmrofel 26/11/2019
  materialsList_.resize(solidRegionNames_.size());
  meshesList_.resize(solidRegionNames_.size());
  dynamicMeshesList_.resize(solidRegionNames_.size());

  forAll(solidRegionNames_, regionI) {
    Info << "=== BEGINNING MATERIAL: " << solidRegionNames_[regionI] << " ===" << endl;
    Info << "Create solid mesh for region " << solidRegionNames_[regionI]
         << " for time = " << runTime.timeName() << nl << endl;

    word& name = solidRegionNames_[regionI];

    IOdictionary matDict_
    (
        IOobject
        (
            solidRegionNames_[regionI]+"/"+solidRegionNames_[regionI]+"Properties",
            runTime.constant(),
            runTime.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Switch debug_(matDict_.lookupOrDefault<Switch>("debug", false));
    Switch movingMesh_="no";
    movingMesh_ = matDict_.lookupOrDefault<Switch>("movingMesh","no");
    if (debug_) {
      Info << "--- movingMesh=" << movingMesh_ << " --- Foam::basicMaterialsModel::basicMaterialsModel"<< endl;
    }
    if(movingMesh_) {
      dynamicMeshesList_.set(
          regionI,
          dynamicFvMesh::New
          (
              IOobject
              (
                  solidRegionNames_[regionI],
                  runTime.timeName(),
                  runTime,
                  IOobject::MUST_READ
              )
          )
      );

      materialsList_.set
      (
          regionI,
          new simpleMaterialModel
          (
              dynamicMeshesList_[regionI],
              solidRegionNames_[regionI]
          )
//          new autoPtr<simpleMaterialModel>
//          (
//              simpleMaterialModel::New
//              (
//                  dynamicMeshesList_[regionI],
//                  word::null,
//                  solidRegionNames_[regionI]
//              )
//          )
      );
    } else {
      meshesList_.set(regionI,
                      new fvMesh
                      (
                          IOobject
                          (
                              solidRegionNames_[regionI],
                              runTime.timeName(),
                              runTime,
                              IOobject::MUST_READ
                          )
                      )
                     );


      materialsList_.set
      (
          regionI,
          new simpleMaterialModel
          (
              meshesList_[regionI],
              solidRegionNames_[regionI]
          )
//          new autoPtr<simpleMaterialModel>
//          (
//              simpleMaterialModel::New
//              (
//                  meshesList_[regionI],
//                  word::null,
//                  solidRegionNames_[regionI]
//              )
//          )
      );
    }
    Info << "=== END MATERIAL: " << solidRegionNames_[regionI] << " ===" << endl << endl;

  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMaterialsModel::~simpleMaterialsModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

void Foam::simpleMaterialsModel::update()
{

  forAll(materialsList_, materialI) {
    materialsList_[materialI].update();  
    simpleMinDeltaTw_ = materialsList_[materialI].updateMinDeltaTw_; // NAmrofel 26/11/2019
  }

  
}
// ************************************************************************* //

