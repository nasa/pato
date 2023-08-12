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
    along with OpenFOAM.  If Fouriert, see <http://www.gnu.org/licenses/>.

Class
    Foam::FourierMaterialPropertiesModel

Description

SourceFiles
    FourierMaterialPropertiesModel.H
    FourierMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef FourierMaterialPropertiesModel_H
#define FourierMaterialPropertiesModel_H

#include "simpleMaterialPropertiesModel.H"
#include "CommonMaterialPropertiesModel.H"
#include "simpleEnergyModel.H"
#include "PMap.H"
#include "MapField.H"

#define bindFun(fun,name) std::bind(&Foam::FourierMaterialPropertiesModel::fun, this, name)
#define newMapField(name,type,init,update) new MapField(name,type,bindFun(init,name),bindFun(update,name))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FourierMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class FourierMaterialPropertiesModel
:
public simpleMaterialPropertiesModel,
public CommonMaterialPropertiesModel
{

 protected:

  //- Number of coefficients
  const scalar& nCoeffs_;

  //- Coefficients
  PMap<scalarList> coeffs_;

  //- Energy Model
  simpleEnergyModel& energyModel_;

  //- Temperature [K]
  volScalarField& T_;

  //- Fourier properties dictionary
  IOdictionary fourierPropDict_;

 public:

  //- Runtime type information
  TypeName("Fourier");

  // Constructors

  //- Construct from mesh and phase name
  FourierMaterialPropertiesModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~FourierMaterialPropertiesModel();

  //- Member Functions

  //- Initialize
  virtual void initialize();

  //- Initialize coeffs
  void initCoeffs(word name);

  //- Initialize k coeffs
  void init_k_coeffs(word name);

  //- Update
  virtual void update();

  //- Update scalar field
  void updateScalarField(word name);

  //- Update tensor field
  void updateTensorField(word name);

  //- Update k
  void update_k(word name);

  //- Fourier scalar function
  void fourierScalarFunction(volScalarField& field, scalarList& coeffs, volScalarField& T);

  //- Fourier tensor function
  void fourierTensorFunction(volTensorField& field, scalarList& coeffs, volScalarField& T);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
