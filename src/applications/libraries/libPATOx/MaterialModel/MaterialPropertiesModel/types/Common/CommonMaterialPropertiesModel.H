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
    along with OpenFOAM.  If Porous_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::CommonMaterialPropertiesModel

Description

SourceFiles
    CommonMaterialPropertiesModel.H
    CommonMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef CommonMaterialPropertiesModel_H
#define CommonMaterialPropertiesModel_H

#include "simpleModel.H"
#include "IOFunctions.H"
#include "PList.H"
#include "PMap.H"
#include "MapField.H"

#define bindFun(fun,name) std::bind(&Foam::CommonMaterialPropertiesModel::fun, this, name)
#define newMapField(name,type,init,update) new MapField(name,type,bindFun(init,name),bindFun(update,name))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CommonMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class CommonMaterialPropertiesModel
{

 protected:

  // Protected Member Variables

  //- simple model
  simpleModel& simpleModel_;

  //- Model name
  const word modelName_;

  //- Type name
  const word typeName_;

  //- List of Map Field object (name,type,init(),update())
  PList<MapField> listMapFields;

  //- Map of Map Field object
  std::map<word, MapField*> mapFields;

  //- Indexes of the fields that will be updated
  List<int> activeFieldIndexes_;

  //- volScalarFields
  PMap<volScalarField> scalarFields_;

  //- volVectorFields
  PMap<volVectorField> vectorFields_;

  //- volTensorFields
  PMap<volTensorField> tensorFields_;

  //- volSymmTensorFields
  PMap<volSymmTensorField> symmTensorFields_;

  //- Switch properties
  PMap<Switch> switchMatProps_;

  //- scalar properties
  PMap<scalar> scalarMatProps_;

  //- vector properties
  PMap<vector> vectorMatProps_;

  //- tensor properties
  PMap<tensor> tensorMatProps_;

  //- dimScalar properties
  PMap<dimensionedScalar> dimScalarMatProps_;

  //- dimVector properties
  PMap<dimensionedVector> dimVectorMatProps_;

  //- dimTensor properties
  PMap<dimensionedTensor> dimTensorMatProps_;

 public:

  // Constructors

  //- Construct from mesh and phase name
  CommonMaterialPropertiesModel(simpleModel& simpleModel, const word& modelName, const word& typeName);

  //- Destructor
  virtual ~CommonMaterialPropertiesModel();

  // Member Functions

  // Initialization functions

  //- Init the porous material properties
  void initialize();

  // Update functions

  //- Update the porous material properties
  void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
