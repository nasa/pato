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

Class
    Foam::simpleModel

Description

SourceFiles
    simpleModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleModel_H
#define simpleModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"
#include <map>
#include "PList.H"
#include "PMap.H"

#define ExitError(msg) FatalError << msg << exit(FatalError)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleModel Declaration
\*---------------------------------------------------------------------------*/

class simpleModel
:
public IOdictionary
{

 private:

  // Private Member Functions

  //- Reference to a word property from this model
  inline word& refWordProp(fileName name);

  //- Reference to a fileName property from this model
  inline fileName& refFileNameProp(fileName name);

  //- Reference to a switch property from this model
  inline Switch& refSwitchProp(fileName name);

  //- Reference to a scalar property from this model
  inline scalar& refScalarProp(fileName name);

  //- Reference to a vector property from this model
  inline vector& refVectorProp(fileName name);

  //- Reference to a tensor property from this model
  inline tensor& refTensorProp(fileName name);

  //- Reference to a dimensionedScalar property from this model
  inline dimensionedScalar& refDimScalarProp(fileName name);

  //- Reference to a dimensionedVector property from this model
  inline dimensionedVector& refDimVectorProp(fileName name);

  //- Reference to a dimensionedTensor property from this model
  inline dimensionedTensor& refDimTensorProp(fileName name);

  //- Replace $PATO_DIR in path
  inline fileName replace_PATO_DIR(fileName name);

 protected:

  // Protected Member Variables

  //- Mesh database
  const fvMesh& mesh_;

  //- Region name
  const word& regionName_;

  //- Model name
  const word& simpleModelName_;

  //- Material dictionary
  IOdictionary materialDict_;

  //- Path of the subdictionary in materialDict_
  const fileName subDictPath_;

  //- Debug Switch
  const Switch debug_;

  //- dynamic mesh flag
  const bool dynamicMesh_;

  //- Material Properties directory
  const fileName materialDirectory_;

  //- constantProperties file name
  const fileName constantPropertiesFile_;

  //- Dictionary constantProperties
  IOdictionary constantPropertiesDictionary_;

  //- Number of mesh patches in this region
  const int nPatches_;

  //- Read fields
  const List<Tuple2<word,word>> createFields_;

  //- Created vol fields names
  wordList createdVolFields_;

  //- Created surface fields names
  wordList createdSurfaceFields_;

  //- Created uniform fields names
  wordList createdUniformFields_;

  //- Switch: model initialized
  Switch modelInitialized_;

  // Static Protected Member Variables

  //- Created word properties names
  static PMap<wordList> wordPropNames_;

  //- Created fileName properties names
  static PMap<wordList> fileNamePropNames_;

  //- Created switch properties names
  static PMap<wordList> switchPropNames_;

  //- Created scalar properties names
  static PMap<wordList> scalarPropNames_;

  //- Created vector properties names
  static PMap<wordList> vectorPropNames_;

  //- Created tensor properties names
  static PMap<wordList> tensorPropNames_;

  //- Created dimensioned scalar properties names
  static PMap<wordList> dimScalarPropNames_;

  //- Created dimensioned vector properties names
  static PMap<wordList> dimVectorPropNames_;

  //- Created dimensioned tensor properties names
  static PMap<wordList> dimTensorPropNames_;

  //- word properties
  static PMap<PList<word>> wordProps_;

  //- fileName properties
  static PMap<PList<fileName>> fileNameProps_;

  //- switch properties
  static PMap<PList<Switch>> switchProps_;

  //- scalar properties
  static PMap<PList<scalar>> scalarProps_;

  //- vector properties
  static PMap<PList<vector>> vectorProps_;

  //- tensor properties
  static PMap<PList<tensor>> tensorProps_;

  //- dimensioned scalar properties
  static PMap<PList<dimensionedScalar>> dimScalarProps_;

  //- dimensioned vector properties
  static PMap<PList<dimensionedVector>> dimVectorProps_;

  //- dimensioned tensor properties
  static PMap<PList<dimensionedTensor>> dimTensorProps_;

 public:

  // Constructors

  //- Construct from mesh and phase name
  simpleModel
  (
      const fvMesh& mesh,
      const word& regionName,
      const word& simpleModelName
  );

  //- Destructor
  virtual ~simpleModel();

  // Public Static Member Variables

  //- Tab Level for Info
  static int tabLevel_;

  // Public Member Functions

  //- Loop over all cells and update all thermodynamic properties
  virtual void update() = 0;

  //- Initialize a word property
  inline word& createWordProp(fileName name, Switch defSwitch="no", word defValue="");

  //- Initialize a fileName property
  inline fileName& createFileNameProp(fileName name, Switch defSwitch="no", fileName defValue="");

  //- Initialize a switch property
  inline Switch& createSwitchProp(fileName name, Switch defSwitch="no", Switch defValue="no");

  //- Initialize a scalar property
  inline scalar& createScalarProp(fileName name, Switch defSwitch="no", scalar defValue=scalar(0));

  //- Initialize a vector property
  inline vector& createVectorProp(fileName name, Switch defSwitch="no", vector defValue=vector::zero);

  //- Initialize a tensor property
  inline tensor& createTensorProp(fileName name, Switch defSwitch="no", tensor defValue=tensor::zero);

  //- Initialize a scalar property
  inline dimensionedScalar& createDimScalarProp
  (fileName name, Switch defSwitch="no", dimensionedScalar defValue=dimensionedScalar("0",dimless,scalar(0)));

  //- Initialize a vector property
  inline dimensionedVector& createDimVectorProp
  (fileName name, Switch defSwitch="no", dimensionedVector defValue=dimensionedVector("0",dimless,vector::zero));

  //- Initialize a tensor property
  inline dimensionedTensor& createDimTensorProp
  (fileName name, Switch defSwitch="no", dimensionedTensor defValue=dimensionedTensor("0",dimless,tensor::zero));

  //- Initialize a vol field (MUST_READ)
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolField(fileName name);

  //- Initialize a surface field (MUST_READ)
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceField(fileName name);

  //- Initialize a uniform field (MUST_READ)
  template<typename Type>
  inline UniformDimensionedField<Type>& createUniformField(fileName name);

  //- Initialize a vol field (READ_IF_PRESENT) from dimensioned value
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolField
  (fileName name, dimensioned<Type> dimValue);

  //- Initialize a surface field (READ_IF_PRESENT) from dimensioned value
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceField
  (fileName name, dimensioned<Type> dimValue);

  //- Initialize a uniform field (READ_IF_PRESENT) from dimensioned value
  template<typename Type>
  inline UniformDimensionedField<Type>& createUniformField
  (fileName name, dimensioned<Type> dimValue);

  //- Initialize a vol field (READ_IF_PRESENT) from dimensioned value
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolField
  (fileName name, dimensioned<Type> dimValue, wordList BC);

  //- Initialize a surface field (READ_IF_PRESENT) from dimensioned value
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceField
  (fileName name, dimensioned<Type> dimValue, wordList BC);

  //- Initialize a vol field (READ_IF_PRESENT) from field
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolField
  (fileName name, const GeometricField<Type, fvPatchField, volMesh>& field);

  //- Initialize a surface field (READ_IF_PRESENT) from field
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceField
  (fileName name, const GeometricField<Type, fvsPatchField, surfaceMesh>& field);

  //- Initialize a vol field (MUST_READ) if not found
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, Switch addField=false);

  //- Initialize a surface field (MUST_READ) if not found
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, Switch addField=false);

  //- Initialize a uniform field (MUST_READ) if not found
  template<typename Type>
  inline UniformDimensionedField<Type>& createUniformFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, Switch addField=false);

  //- Initialize a vol field (READ_IF_PRESENT) from dimensioned value if not found
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, Switch addField=false);

  //- Initialize a surface field (READ_IF_PRESENT) from dimensioned value if not found
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, Switch addField=false);

  //- Initialize a uniform field (READ_IF_PRESENT) from dimensioned value if not found
  template<typename Type>
  inline UniformDimensionedField<Type>& createUniformFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, Switch addField=false);

  //- Initialize a vol field (READ_IF_PRESENT) from dimensioned value if not found
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, wordList BC, Switch addField=false);

  //- Initialize a surface field (READ_IF_PRESENT) from dimensioned value if not found
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, wordList BC, Switch addField=false);

  //- Initialize a vol field (READ_IF_PRESENT) from field if not found
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& createVolFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, const GeometricField<Type, fvPatchField, volMesh>& field, Switch addField=false);

  //- Initialize a surface field (READ_IF_PRESENT) from field if not found
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& createSurfaceFieldIfNotFound
  (const simpleModel& constOtherModel_, fileName name, const GeometricField<Type, fvsPatchField, surfaceMesh>& field, Switch addField=false);

  //- Increase the tab level. To be called at the beginning of the model initinalization after parent class initialization.
  inline Switch startModelInit();

  //- Decrease the tab level. To be called at the end of the model initinalization.
  inline void modelInitialized();

  // Public Member Constant Functions

  //- Mesh database
  inline const fvMesh& mesh() const;

  //- Debug flag
  inline const Switch& debug() const;

  //- Region name
  inline const word& regionName() const;

  //- Material dictionary
  inline const IOdictionary& materialDict() const;

  //- Model name
  inline const word& simpleModelName() const;

  //- Created vol fields
  inline const wordList& createdVolFields() const;

  //- Created surface fields
  inline const wordList& createdSurfaceFields() const;

  //- Created uniform fields
  inline const wordList& createdUniformFields() const;

  //- Reference to another model
  template<typename Type>
  inline Type& refModel() const;

  //- Reference to a PList of vol field from this model
  template<typename Type>
  inline PList<GeometricField<Type, fvPatchField, volMesh>> refVolFieldPList(fileName name, int startIndexName=0, wordList names=wordList(), Switch acceptEmpty=false) const;

  //- Reference to a vol field from this model
  template<typename Type>
  inline GeometricField<Type, fvPatchField, volMesh>& refVolField(fileName name) const;

  //- Reference to a vol/surface field from this model
  template<typename Type>
  inline GeometricField<Type, fvsPatchField, surfaceMesh>& refSurfaceField(fileName name) const;

  //- Reference to a vol/surface field from this model
  template<typename Type>
  inline UniformDimensionedField<Type>& refUniformField(fileName name) const;

  //- Return true if found the volField name in this model
  inline Switch foundVolField(fileName name) const;

  //- Return true if found the surfaceField name in this model
  inline Switch foundSurfaceField(fileName name) const;

  //- Return true if found the uniformField name in this model
  inline Switch foundUniformField(fileName name) const;

  // Static Public Member Functions

  //- Get the tab level
  inline static const word getTabLevel(word end=" ");

};

#include "simpleModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
