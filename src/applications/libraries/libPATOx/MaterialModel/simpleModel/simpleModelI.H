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

#include "simpleModel.H"

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

inline Foam::word& Foam::simpleModel::refWordProp(Foam::fileName name)
{
  int index = -1;
  forAll(wordPropNames_[regionName_], elemI) {
    if(wordPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " wordProp" << endl;
  return wordProps_[regionName_][index];
}

inline Foam::fileName& Foam::simpleModel::refFileNameProp(Foam::fileName name)
{
  int index = -1;
  forAll(fileNamePropNames_[regionName_], elemI) {
    if(fileNamePropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " wordProp" << endl;
  return fileNameProps_[regionName_][index];
}

inline Foam::Switch& Foam::simpleModel::refSwitchProp(Foam::fileName name)
{
  int index = -1;
  forAll(switchPropNames_[regionName_], elemI) {
    if(switchPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " switchProp" << endl;
  return switchProps_[regionName_][index];
}

inline Foam::scalar& Foam::simpleModel::refScalarProp(Foam::fileName name)
{
  int index = -1;
  forAll(scalarPropNames_[regionName_], elemI) {
    if(scalarPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " scalarProp" << endl;
  return scalarProps_[regionName_][index];
}

inline Foam::vector& Foam::simpleModel::refVectorProp(Foam::fileName name)
{
  int index = -1;
  forAll(vectorPropNames_[regionName_], elemI) {
    if(vectorPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " vectorProp" << endl;
  return vectorProps_[regionName_][index];
}

inline Foam::tensor& Foam::simpleModel::refTensorProp(Foam::fileName name)
{
  int index = -1;
  forAll(tensorPropNames_[regionName_], elemI) {
    if(tensorPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " tensorProp" << endl;
  return tensorProps_[regionName_][index];
}

inline Foam::dimensionedScalar& Foam::simpleModel::refDimScalarProp(Foam::fileName name)
{
  int index = -1;
  forAll(dimScalarPropNames_[regionName_], elemI) {
    if(dimScalarPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " dimScalarProp" << endl;
  return dimScalarProps_[regionName_][index];
}

inline Foam::dimensionedVector& Foam::simpleModel::refDimVectorProp(Foam::fileName name)
{
  int index = -1;
  forAll(dimVectorPropNames_[regionName_], elemI) {
    if(dimVectorPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " dimVectorProp" << endl;
  return dimVectorProps_[regionName_][index];
}

inline Foam::dimensionedTensor& Foam::simpleModel::refDimTensorProp(Foam::fileName name)
{
  int index = -1;
  forAll(dimTensorPropNames_[regionName_], elemI) {
    if(dimTensorPropNames_[regionName_][elemI]==name) {
      index = elemI;
      break;
    }
  }
  if (index<0) {
    FatalError << name << " property not found in region " << regionName_ << " using "
               << simpleModelName_ << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " dimTensorProp" << endl;
  return dimTensorProps_[regionName_][index];
}

inline Foam::fileName Foam::simpleModel::replace_PATO_DIR(fileName name)
{
  fileName new_path = name.replaceAll(getEnv("PATO_DIR"),"$PATO_DIR");
  new_path = new_path.replaceAll(getEnv("PATO_TUTORIALS"),"$PATO_TUTORIALS");
  return new_path;
}

// * * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * * //

inline Foam::word& Foam::simpleModel::createWordProp(Foam::fileName name, Foam::Switch defSwitch, Foam::word defValue)
{
  if (foundInList((word)name,wordPropNames_[regionName_])) {
    return refWordProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    wordProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(subDictPath_) << ": ";
      wordProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " scalar property with default value: ";
        wordProps_[regionName_].append(new word(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  wordPropNames_[regionName_].append(name);
  const int wordProps_size = wordProps_[regionName_].size();
  Info << (word) name << "=" << wordProps_[regionName_][wordProps_size-1] << endl;
  return wordProps_[regionName_][wordProps_size-1];
}

inline Foam::fileName& Foam::simpleModel::createFileNameProp(Foam::fileName name, Foam::Switch defSwitch, Foam::fileName defValue)
{
  if (foundInList((word)name,fileNamePropNames_[regionName_])) {
    return refFileNameProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    fileNameProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(subDictPath_) << ": ";
      fileNameProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " scalar property with default value: ";
        fileNameProps_[regionName_].append(new fileName(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  fileNamePropNames_[regionName_].append(name);
  const int fileNameProps_size = fileNameProps_[regionName_].size();
  Info << (word) name << "=" << fileNameProps_[regionName_][fileNameProps_size-1] << endl;
  return fileNameProps_[regionName_][fileNameProps_size-1];
}

inline Foam::Switch& Foam::simpleModel::createSwitchProp(Foam::fileName name, Foam::Switch defSwitch, Foam::Switch defValue)
{
  if (foundInList((word)name,switchPropNames_[regionName_])) {
    return refSwitchProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    switchProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(subDictPath_) << ": ";
      switchProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " scalar property with default value: ";
        switchProps_[regionName_].append(new Switch(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  switchPropNames_[regionName_].append(name);
  const int switchProps_size = switchProps_[regionName_].size();
  Info << (word) name << "=" << switchProps_[regionName_][switchProps_size-1] << endl;
  return switchProps_[regionName_][switchProps_size-1];
}

inline Foam::scalar& Foam::simpleModel::createScalarProp(Foam::fileName name, Foam::Switch defSwitch, Foam::scalar defValue)
{
  if (foundInList((word)name,scalarPropNames_[regionName_])) {
    return refScalarProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    scalarProps_[regionName_].append(new scalar(readScalar(constantPropertiesDictionary_.lookup(name))));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " scalar property from " << replace_PATO_DIR(subDictPath_) << ": ";
      scalarProps_[regionName_].append(new scalar(readScalar(materialDict_.subDict(simpleModelName_).lookup(name))));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " scalar property with default value: ";
        scalarProps_[regionName_].append(new scalar(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  scalarPropNames_[regionName_].append(name);
  const int scalarProps_size = scalarProps_[regionName_].size();
  Info << (word) name << "=" << scalarProps_[regionName_][scalarProps_size-1] << endl;
  return scalarProps_[regionName_][scalarProps_size-1];
}

inline Foam::vector& Foam::simpleModel::createVectorProp(Foam::fileName name, Foam::Switch defSwitch, Foam::vector defValue)
{
  if (foundInList((word)name,vectorPropNames_[regionName_])) {
    return refVectorProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " vector property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    vectorProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " vector property from " << replace_PATO_DIR(subDictPath_) << ": ";
      vectorProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " vector property with default value: ";
        vectorProps_[regionName_].append(new vector(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  vectorPropNames_[regionName_].append(name);
  const int vectorProps_size=vectorProps_[regionName_].size();
  Info << (word) name << "=" << vectorProps_[regionName_][vectorProps_size-1] << endl;
  return vectorProps_[regionName_][vectorProps_size-1];
}

inline Foam::tensor& Foam::simpleModel::createTensorProp(Foam::fileName name, Foam::Switch defSwitch, Foam::tensor defValue)
{
  if (foundInList((word)name,tensorPropNames_[regionName_])) {
    return refTensorProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " tensor property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    tensorProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " tensor property from " << replace_PATO_DIR(subDictPath_) << ": ";
      tensorProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " tensor property with default value: ";
        tensorProps_[regionName_].append(new tensor(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  tensorPropNames_[regionName_].append(name);
  const int tensorProps_size=tensorProps_[regionName_].size();
  Info << (word) name << "=" << tensorProps_[regionName_][tensorProps_size-1] << endl;
  return tensorProps_[regionName_][tensorProps_size-1];
}

inline Foam::dimensionedScalar& Foam::simpleModel::createDimScalarProp(Foam::fileName name, Foam::Switch defSwitch, Foam::dimensionedScalar defValue)
{
  if (foundInList((word)name,dimScalarPropNames_[regionName_])) {
    return refDimScalarProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " dimScalar property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    dimScalarProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " dimScalar property from " << replace_PATO_DIR(subDictPath_) << ": ";
      dimScalarProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " dimScalar property with default value: ";
        dimScalarProps_[regionName_].append(new dimensionedScalar(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  dimScalarPropNames_[regionName_].append(name);
  const int dimScalarProps_size=dimScalarProps_[regionName_].size();
  Info << (word) name << "=" << dimScalarProps_[regionName_][dimScalarProps_size-1] << endl;
  return dimScalarProps_[regionName_][dimScalarProps_size-1];
}

inline Foam::dimensionedVector& Foam::simpleModel::createDimVectorProp(Foam::fileName name, Foam::Switch defSwitch, Foam::dimensionedVector defValue)
{
  if (foundInList((word)name,dimVectorPropNames_[regionName_])) {
    return refDimVectorProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " dimVector property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    dimVectorProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " dimVector property from " << replace_PATO_DIR(subDictPath_) << ": ";
      dimVectorProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " dimVector property with default value: ";
        dimVectorProps_[regionName_].append(new dimensionedVector(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  dimVectorPropNames_[regionName_].append(name);
  const int dimVectorProps_size=dimVectorProps_[regionName_].size();
  Info << (word) name << "=" << dimVectorProps_[regionName_][dimVectorProps_size-1] << endl;
  return dimVectorProps_[regionName_][dimVectorProps_size-1];
}

inline Foam::dimensionedTensor& Foam::simpleModel::createDimTensorProp(Foam::fileName name, Foam::Switch defSwitch, Foam::dimensionedTensor defValue)
{
  if (foundInList((word)name,dimTensorPropNames_[regionName_])) {
    return refDimTensorProp(name);
  }
  if (constantPropertiesDictionary_.found(name)) {
    Info << getTabLevel() << "Reading " << name << " dimTensor property from " << replace_PATO_DIR(constantPropertiesFile_) << ": ";
    dimTensorProps_[regionName_].append(constantPropertiesDictionary_.lookup(name));
  } else {
    if (materialDict_.subDict(simpleModelName_).found(name)) {
      Info << getTabLevel() << "Reading " << name << " dimTensor property from " << replace_PATO_DIR(subDictPath_) << ": ";
      dimTensorProps_[regionName_].append(materialDict_.subDict(simpleModelName_).lookup(name));
    } else {
      if (defSwitch) {
        Info << getTabLevel() << "Create " << name << " dimTensor property with default value: ";
        dimTensorProps_[regionName_].append(new dimensionedTensor(defValue));
      } else {
        FatalError << simpleModelName_ << "Model: "<< (fileName) name
                   << " property not found in " << constantPropertiesFile_
                   << " or " << subDictPath_
                   << exit(FatalError);
      }
    }
  }
  dimTensorPropNames_[regionName_].append(name);
  const int dimTensorProps_size=dimTensorProps_[regionName_].size();
  Info << (word) name << "=" << dimTensorProps_[regionName_][dimTensorProps_size-1] << endl;
  return dimTensorProps_[regionName_][dimTensorProps_size-1];
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolField(Foam::fileName name)
{
  Info << getTabLevel() << "Create " << name << " volField (MUST_READ)" << endl;
  createdVolFields_.append(name);
  return meshUniqueConstruct<Type, fvPatchField, volMesh>(mesh_,name,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceField(Foam::fileName name)
{
  Info << getTabLevel() << "Create " << name << " surfaceField (MUST_READ)" << endl;
  createdSurfaceFields_.append(name);
  return meshUniqueConstruct<Type, fvsPatchField, surfaceMesh>(mesh_,name,simpleModelName_);
}

template<typename Type>
inline UniformDimensionedField<Type>& Foam::simpleModel::createUniformField(Foam::fileName name)
{
  Info << getTabLevel() << "Create " << name << " uniformField (MUST_READ)" << endl;
  createdUniformFields_.append(name);
  return meshUniqueConstruct<Type>(mesh_,name,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolField(Foam::fileName name, Foam::dimensioned<Type> dimValue)
{
  Info << getTabLevel() << "Create " << name << " volField (READ_IF_PRESENT)" << endl;
  createdVolFields_.append(name);
  return meshUniqueConstruct<Type, fvPatchField, volMesh>(mesh_,name,dimValue,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceField(Foam::fileName name, Foam::dimensioned<Type> dimValue)
{
  Info << getTabLevel() << "Create " << name << " surfaceField (READ_IF_PRESENT)" << endl;
  createdSurfaceFields_.append(name);
  return meshUniqueConstruct<Type, fvsPatchField, surfaceMesh>(mesh_,name,dimValue,simpleModelName_);
}

template<typename Type>
inline UniformDimensionedField<Type>& Foam::simpleModel::createUniformField(Foam::fileName name, Foam::dimensioned<Type> dimValue)
{
  Info << getTabLevel() << "Create " << name << " uniformField (READ_IF_PRESENT)" << endl;
  createdUniformFields_.append(name);
  return meshUniqueConstruct<Type>(mesh_,name,dimValue,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolField(Foam::fileName name, Foam::dimensioned<Type> dimValue, Foam::wordList BC)
{
  Info << getTabLevel() << "Create " << name << " volField (READ_IF_PRESENT)" << endl;
  createdVolFields_.append(name);
  return meshUniqueConstruct<Type, fvPatchField, volMesh>(mesh_,name,dimValue,BC,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceField(Foam::fileName name, Foam::dimensioned<Type> dimValue, Foam::wordList BC)
{
  Info << getTabLevel() << "Create " << name << " surfaceField (READ_IF_PRESENT)" << endl;
  createdSurfaceFields_.append(name);
  return meshUniqueConstruct<Type, fvsPatchField, surfaceMesh>(mesh_,name,dimValue,BC,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolField(Foam::fileName name, const GeometricField<Type, fvPatchField, volMesh>& field)
{
  createdVolFields_.append(name);
  return meshUniqueConstruct<Type, fvPatchField, volMesh>(mesh_,name,field,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceField(Foam::fileName name, const GeometricField<Type, fvsPatchField, surfaceMesh>& field)
{
  createdSurfaceFields_.append(name);
  return meshUniqueConstruct<Type, fvsPatchField, surfaceMesh>(mesh_,name,field,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundVolField(name)) {
    if (addField) createdVolFields_.append(name);
    return otherModel.refVolField<Type>(name);
  } else {
    return createVolField<Type>(name);
  }
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundSurfaceField(name)) {
    if (addField) createdSurfaceFields_.append(name);
    return otherModel.refSurfaceField<Type>(name);
  } else {
    return createSurfaceField<Type>(name);
  }
}

template<typename Type>
inline UniformDimensionedField<Type>& Foam::simpleModel::createUniformFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundUniformField(name)) {
    if (addField) createdUniformFields_.append(name);
    return otherModel.refUniformField<Type>(name);
  } else {
    return createUniformField<Type>(name);
  }
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundVolField(name)) {
    if (addField) createdVolFields_.append(name);
    return otherModel.refVolField<Type>(name);
  } else {
    return createVolField<Type>(name,dimValue);
  }
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundSurfaceField(name)) {
    if (addField) createdSurfaceFields_.append(name);
    return otherModel.refSurfaceField<Type>(name);
  } else {
    return createSurfaceField<Type>(name,dimValue);
  }
}

template<typename Type>
inline UniformDimensionedField<Type>& Foam::simpleModel::createUniformFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundUniformField(name)) {
    if (addField) createdUniformFields_.append(name);
    return otherModel.refUniformField<Type>(name);
  } else {
    return createUniformField<Type>(name,dimValue);
  }
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, wordList BC, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundVolField(name)) {
    if (addField) createdVolFields_.append(name);
    return otherModel.refVolField<Type>(name);
  } else {
    return createVolField<Type>(name,dimValue,BC);
  }
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceFieldIfNotFound(const simpleModel& constOtherModel_, fileName name, dimensioned<Type> dimValue, wordList BC, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundSurfaceField(name)) {
    if (addField) createdSurfaceFields_.append(name);
    return otherModel.refSurfaceField<Type>(name);
  } else {
    return createSurfaceField<Type>(name,dimValue,BC);
  }
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::createVolFieldIfNotFound
(const simpleModel& constOtherModel_, fileName name, const GeometricField<Type, fvPatchField, volMesh>& field, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundVolField(name)) {
    if (addField) createdVolFields_.append(name);
    return otherModel.refVolField<Type>(name);
  } else {
    return createVolField<Type>(name,field);
  }
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::createSurfaceFieldIfNotFound
(const simpleModel& constOtherModel_, fileName name, const GeometricField<Type, fvsPatchField, surfaceMesh>& field, Switch addField)
{
  simpleModel& otherModel=const_cast<simpleModel&>(constOtherModel_);
  if (otherModel.foundSurfaceField(name)) {
    if (addField) createdSurfaceFields_.append(name);
    return otherModel.refVolField<Type>(name);
  } else {
    return createSurfaceField<Type>(name,field);
  }
}

inline Switch Foam::simpleModel::startModelInit()
{
  tabLevel_++;
  modelInitialized_=false;
  return modelInitialized_;
}

inline void Foam::simpleModel::modelInitialized()
{
  if (!modelInitialized_) tabLevel_--;
  modelInitialized_="yes";
}

// * * * * * * * * * * * * * * * * Public Member Constant Functions  * * * * * * * * * * * * * * //

template<typename Type>
inline Type& Foam::simpleModel::refModel() const
{
  return meshLookupOrConstructModel<Type>(mesh_,regionName_,Type::modelName);
}

template<typename Type>
inline PList<GeometricField<Type, fvPatchField, volMesh>> Foam::simpleModel::refVolFieldPList(Foam::fileName name, int startIndexName, wordList names, Switch acceptEmpty) const
{
  const int size_names=names.size();
  int index=startIndexName;
  word filedName=name;
  PList<GeometricField<Type, fvPatchField, volMesh>> list;
  while(true) {
    string fieldName=name+"["+Foam::name(index)+"]";
    if (size_names!=0) {
      if (index >= size_names) {
        break;
      } else {
        fieldName=name+"["+names[index]+"]";
      }
    }
    if (!foundInList((word)fieldName,createdVolFields_)) {
      break;
    }
    Info << getTabLevel() << "Reference to " << fieldName << " volField from " << simpleModelName_ << endl;
    list.append(meshLookup<Type, fvPatchField, volMesh>(mesh_,fieldName,simpleModelName_));
    index++;
  }
  if (list.size()==0) {
    if (acceptEmpty) {
      Info << getTabLevel() << "Reference to " << name << " empty PList<volField> from " << simpleModelName_ << endl;
    } else {
      FatalError << name << "[0] field not found in region "
      << regionName_ << " of " << simpleModelName_
      << "Model" << exit(FatalError);
    }
  }
  return list;
}

template<typename Type>
inline GeometricField<Type, fvPatchField, volMesh>& Foam::simpleModel::refVolField(Foam::fileName name) const
{
  if (!foundInList((word)name,createdVolFields_)) {
    FatalError << name << " field not found in region "
               << regionName_ << " of " << simpleModelName_
               << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " volField from " << simpleModelName_ << endl;
  return meshLookup<Type, fvPatchField, volMesh>(mesh_,name,simpleModelName_);
}

template<typename Type>
inline GeometricField<Type, fvsPatchField, surfaceMesh>& Foam::simpleModel::refSurfaceField(Foam::fileName name) const
{
  if (!foundInList((word)name,createdSurfaceFields_)) {
    FatalError << name << " field not found in region "
               << regionName_ << " of " << simpleModelName_
               << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " surfaceField from " << simpleModelName_ << endl;
  return meshLookup<Type, fvsPatchField, surfaceMesh>(mesh_,name,simpleModelName_);
}

template<typename Type>
inline UniformDimensionedField<Type>& Foam::simpleModel::refUniformField(Foam::fileName name) const
{
  if (!foundInList((word)name,createdUniformFields_)) {
    FatalError << name << " field not found in region "
               << regionName_ << " of " << simpleModelName_
               << "Model" << exit(FatalError);
  }
  Info << getTabLevel() << "Reference to " << name << " uniformField from " << simpleModelName_ << endl;
  return meshLookup<UniformDimensionedField<Type>>(mesh_,name,simpleModelName_);
}

inline Switch Foam::simpleModel::foundVolField(Foam::fileName name) const
{
  return foundInList((word)name,createdVolFields_);
}

inline Switch Foam::simpleModel::foundSurfaceField(Foam::fileName name) const
{
  return foundInList((word)name,createdSurfaceFields_);
}

inline Switch Foam::simpleModel::foundUniformField(Foam::fileName name) const
{
  return foundInList((word)name,createdUniformFields_);
}

inline const Foam::fvMesh& Foam::simpleModel::mesh() const
{
  return mesh_;
}

inline const Foam::Switch& Foam::simpleModel::debug() const
{
  return debug_;
}

inline const Foam::word& Foam::simpleModel::regionName() const
{
  return regionName_;
}

inline const Foam::IOdictionary& Foam::simpleModel::materialDict() const
{
  return materialDict_;
}

inline const Foam::word& Foam::simpleModel::simpleModelName() const
{
  return simpleModelName_;
}

inline const Foam::wordList& Foam::simpleModel::createdVolFields() const
{
  return createdVolFields_;
}

inline const Foam::wordList& Foam::simpleModel::createdSurfaceFields() const
{
  return createdSurfaceFields_;
}

inline const Foam::wordList& Foam::simpleModel::createdUniformFields() const
{
  return createdUniformFields_;
}

// * * * * * * * * * * * * * * * * Static Public Member Functions  * * * * * * * * * * * * * * //

inline const Foam::word Foam::simpleModel::getTabLevel(word end)
{
  if (tabLevel_<=0) {
    tabLevel_=0;
    return end;
  }
  word sign="|--";
  const int sign_size=sign.size();
  word tab="";
  for (int i = 0; i < sign_size; i++) {
    if (i==0) {
      tab+=".";
    } else {
      tab+=" ";
    }
  }
  word tot_tab="";
  for(int i = 0; i < tabLevel_; i++) {
    tot_tab+=tab;
  }
  word out=tot_tab+sign+end;
  return out;
}

// ************************************************************************* //

