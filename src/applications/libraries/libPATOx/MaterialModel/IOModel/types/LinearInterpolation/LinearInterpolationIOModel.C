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
    along with OpenFOAM.  If LinearInterpolationt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LinearInterpolationIOModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LinearInterpolationIOModel::LinearInterpolationIOModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
noIOModel(mesh, regionName),
startModelInit_(startModelInit()),
outputList_(simpleIOModel::outputList_),
setLinear_(materialDict_.subDict("IO").lookup("setLinear")),
linearFieldList_(materialDict_.subDict("IO").lookup("linearFieldList")),
initOutput_(simpleIOModel::initOutput())
{
  if (setLinear_) {
    forAll(linearFieldList_, fieldI) {
      if (!mesh_.objectRegistry::foundObject<volScalarField>(linearFieldList_[fieldI])) {
        FatalErrorInFunction << linearFieldList_[fieldI] << " field not found in mesh. Linear interpolation problem in Output model." << exit(FatalError);
      }
      volScalarField& field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(linearFieldList_[fieldI]));
      scalar field_top(readScalar(simpleIOModel::materialDict_.subDict("IO").lookup(linearFieldList_[fieldI]+"_top")));
      scalar field_interface(readScalar(simpleIOModel::materialDict_.subDict("IO").lookup(linearFieldList_[fieldI]+"_interface")));
      scalar field_bottom(readScalar(simpleIOModel::materialDict_.subDict("IO").lookup(linearFieldList_[fieldI]+"_bottom")));
      word component = simpleIOModel::materialDict_.subDict("IO").lookup("component");
      int indexComponent = -1;
      if (component=="x") {
        indexComponent=0;
      } else if (component == "y") {
        indexComponent=1;
      } else if (component == "z") {
        indexComponent=2;
      } else {
        FatalErrorInFunction << "component has to be: x,y or z." << exit(FatalError);
      }
      scalar y_top(readScalar(simpleIOModel::materialDict_.subDict("IO").lookup(component+"_top")));
      scalar y_interface(readScalar(simpleIOModel::materialDict_.subDict("IO").lookup(component+"_interface")));
      scalar y_bottom(readScalar(simpleIOModel::materialDict_.subDict("IO").lookup(component+"_bottom")));

      // initialize linear field
      forAll(field, celli) {
        if (mesh.C()[celli][indexComponent]>y_interface) {
          field[celli] = field_interface + (field_top-field_interface)*(mesh.C()[celli][indexComponent]-y_interface)/(y_top-y_interface);
        } else {
          field[celli] = field_bottom + (field_interface-field_bottom)*(mesh.C()[celli][indexComponent]-y_bottom)/(y_interface-y_bottom);
        }
      }
      field.write();
    }
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LinearInterpolationIOModel::~LinearInterpolationIOModel()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LinearInterpolationIOModel::update()
{
  noIOModel::writeOutput();
}

// ************************************************************************* //
