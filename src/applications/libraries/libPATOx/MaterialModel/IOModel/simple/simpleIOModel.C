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

#include "simpleIOModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleIOModel, 0);
defineRunTimeSelectionTable(simpleIOModel, fvMesh);
}

/* * * * * * * * * * * * * * * public static members * * * * * * * * * * * * * */

const Foam::word Foam::simpleIOModel::modelName="IO";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleIOModel::simpleIOModel
(
    const fvMesh& mesh,
    const word& regionName
):
simpleModel(mesh,regionName,modelName),
#if defined(FOAM_EXTEND)
foam_extend_mesh_(mesh_.time().db().lookupObject<simpleMaterialsModel>("MaterialsModel").get_foam_extend_mesh(regionName)),
#endif
phaseName_(""),
infoDebug_(materialDict_.isDict("IO")?materialDict_.subDict("IO").lookupOrDefault<label>("infoDebug", lduMatrix::debug): lduMatrix::debug),
materialDictPath_(mesh.time().constant()+"/"+regionName+"/" + IOobject::groupName(regionName+"Properties", phaseName_)),
readFilesList_(materialDict_.isDict("IO")?materialDict_.subDict("IO").lookupOrDefault<List<fileName> >("readFiles",nullListFileName):nullListFileName),
writeFields_(mesh_.time().controlDict().lookupOrDefault<Switch>("writeFields","yes"))
{
  if (infoDebug_ != lduMatrix::debug) {
    // IO of the solvers
    lduMatrix::debug=infoDebug_;
    solverPerformance::debug=infoDebug_;
    Info << "lduMatrix::debug=" << infoDebug_ << endl;
    Info << "solverPerformance::debug=" << infoDebug_ << endl;
  }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleIOModel> Foam::simpleIOModel::New
(
    const fvMesh& mesh,
    const word& regionName
)
{
  IOdictionary MaterialDict
  (
      IOobject
      (
          regionName+"Properties",
          mesh.time().constant(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE,
          false
      )
  );

  word IOModelTypeName= "no";

  if (MaterialDict.isDict("IO")) {
    dictionary IOModelDict = MaterialDict.subDict("IO");

    IOModelTypeName =
        word(IOModelDict.lookupOrDefault<word>("IOType","no"));

  }

  typename simpleIOModel::fvMeshConstructorTable::iterator cstrIter =
      simpleIOModel::fvMeshConstructorTablePtr_->find(IOModelTypeName);

  if (cstrIter == simpleIOModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleIOModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleIOModel::typeName << " type "
        << IOModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleIOModel>(cstrIter()(mesh, regionName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleIOModel::~simpleIOModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

bool Foam::simpleIOModel::initOutput()
{

  if (materialDict_.isDict("IO")) {

    forAll(readFilesList_, fileI) {
      filesData_.append(readFileData(readFilesList_[fileI]));
    }

    wordList outputList_tmp(materialDict_.subDict("IO").lookupOrDefault<wordList>("writeFields",nullList));
    wordList probingDictNames_tmp(materialDict_.subDict("IO").lookupOrDefault<wordList>("probingFunctions",nullList));
    if (!writeFields_) {
      outputList_tmp = wordList();
      const_cast<fvMesh&>(mesh_).objectRegistry::writeOpt()=IOobject::NO_WRITE;
    }
    outputList_=outputList_tmp;
    probingDictNames_=probingDictNames_tmp;

#if defined(FOAM_EXTEND)
    initFoamExtendFields();
#endif

#if defined(darwin64)
    // Mute the libsampling.so warning
    const int Warning_level = Warning.level;
    Warning.level=0;
#endif
    forAll(probingDictNames_, nameI) {
      sampleFunctions_.append(new SampleFunction(mesh_, phaseName_, regionName_, probingDictNames_[nameI]));
    }
#if defined(darwin64)
    Warning.level=Warning_level;
#endif

    foundFieldsInMeshAll(mesh_,outputList_);
  }
  return true;
}

#if defined(FOAM_EXTEND)

#define meshFoundObject(mesh,type,name) mesh.objectRegistry::foundObject<type>(name)
#define meshLookupObject(mesh,type,name) mesh.objectRegistry::lookupObject<type>(name)
#define createDimensions(field) Foam_extend_::dimensionSet ds = field.dimensions();                               \
 dimensionSet foam_dimensionSet(ds[Foam_extend_::dimensionSet::MASS],ds[Foam_extend_::dimensionSet::LENGTH],      \
                                ds[Foam_extend_::dimensionSet::TIME],ds[Foam_extend_::dimensionSet::TEMPERATURE], \
                                ds[Foam_extend_::dimensionSet::MOLES],ds[Foam_extend_::dimensionSet::CURRENT],    \
                                ds[Foam_extend_::dimensionSet::LUMINOUS_INTENSITY])

void Foam::simpleIOModel::initFoamExtendFields()
{
  // Search fields in OpenFOAM mesh database
  wordList fields_not_found;
  forAll(outputList_, wordI) {
    bool found = false;
    word name = outputList_[wordI];
    if (meshFoundObject(mesh_,volScalarField,name) || meshFoundObject(mesh_,volVectorField,name) || \
        meshFoundObject(mesh_,volTensorField,name) || meshFoundObject(mesh_,surfaceScalarField,name) || \
        meshFoundObject(mesh_,surfaceVectorField,name) || meshFoundObject(mesh_,surfaceTensorField,name) ) {
      found = true;
    }
    if (found) continue;
    fields_not_found.append(outputList_[wordI]);
  }

  // Search fields in Foam Extend mesh database
  forAll(fields_not_found, wordI) {
    word name = fields_not_found[wordI];
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volScalarField,name) || meshFoundObject(foam_extend_mesh_,Foam_extend_::volVectorField,name) || \
        meshFoundObject(foam_extend_mesh_,Foam_extend_::volTensorField,name) || meshFoundObject(foam_extend_mesh_,Foam_extend_::volSymmTensorField,name) || \
        meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceScalarField,name) || meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceVectorField,name) || \
        meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceTensorField,name) || meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceSymmTensorField,name)) {
      foam_extend_write_fields_.append(fields_not_found[wordI]);
    }
  }

  // Create the Foam Extend fields in OpenFOAM mesh database
  forAll(foam_extend_write_fields_, wordI) {
    word name = foam_extend_write_fields_[wordI];
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volScalarField,name)) {
      const Foam_extend_::volScalarField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volScalarField,name));
      createDimensions(foam_extend_field);
      volScalarField& field(meshUniqueConstruct<scalar, fvPatchField, volMesh>(mesh_, name.c_str(), dimensionedScalar("0",foam_dimensionSet,0), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volVectorField,name)) {
      const Foam_extend_::volVectorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volVectorField,name));
      createDimensions(foam_extend_field);
      volVectorField& field(meshUniqueConstruct<vector, fvPatchField, volMesh>(mesh_, name.c_str(), \
                            dimensioned<vector>("0",foam_dimensionSet,vector(0,0,0)), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volTensorField,name)) {
      const Foam_extend_::volTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volTensorField,name));
      createDimensions(foam_extend_field);
      volTensorField& field(meshUniqueConstruct<tensor, fvPatchField, volMesh>(mesh_, name.c_str(), \
                            dimensioned<tensor>("0",foam_dimensionSet,tensor(0,0,0,0,0,0,0,0,0)), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volSymmTensorField,name)) {
      const Foam_extend_::volSymmTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volSymmTensorField,name));
      createDimensions(foam_extend_field);
      volSymmTensorField& field(meshUniqueConstruct<symmTensor, fvPatchField, volMesh>(mesh_, name.c_str(), \
                                dimensioned<symmTensor>("0",foam_dimensionSet,symmTensor(0,0,0,0,0,0)), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceScalarField,name)) {
      const Foam_extend_::surfaceScalarField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceScalarField,name));
      createDimensions(foam_extend_field);
      surfaceScalarField& field(meshUniqueConstruct<scalar, fvsPatchField, surfaceMesh>(mesh_, name.c_str(), \
                                dimensionedScalar("0",foam_dimensionSet,0), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceVectorField,name)) {
      const Foam_extend_::surfaceVectorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceVectorField,name));
      createDimensions(foam_extend_field);
      surfaceVectorField& field(meshUniqueConstruct<vector, fvsPatchField, surfaceMesh>(mesh_, name.c_str(), \
                                dimensioned<vector>("0",foam_dimensionSet,vector(0,0,0)), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceTensorField,name)) {
      const Foam_extend_::surfaceTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceTensorField,name));
      createDimensions(foam_extend_field);
      surfaceTensorField& field(meshUniqueConstruct<tensor, fvsPatchField, surfaceMesh>(mesh_, name.c_str(), \
                                dimensioned<tensor>("0",foam_dimensionSet,tensor(0,0,0,0,0,0,0,0,0)), modelName, IOobject::NO_READ));
    }
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceSymmTensorField,name)) {
      const Foam_extend_::surfaceSymmTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceSymmTensorField,name));
      createDimensions(foam_extend_field);
      surfaceSymmTensorField& field(meshUniqueConstruct<symmTensor, fvsPatchField, surfaceMesh>(mesh_, name.c_str(), \
                                    dimensioned<symmTensor>("0",foam_dimensionSet,symmTensor(0,0,0,0,0,0)), modelName, IOobject::NO_READ));
    }
  }
}

void Foam::simpleIOModel::updateFoamExtendFields()
{
  // Update all the Foam Extend write fields
  forAll(foam_extend_write_fields_, wordI) {
    word name = foam_extend_write_fields_[wordI];
    // volScalarField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volScalarField,name)) {
      const Foam_extend_::volScalarField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volScalarField,name));
      volScalarField& field(const_cast<volScalarField&>(meshLookupObject(mesh_,volScalarField,name)));
      forAll(field, cellI) {
        field[cellI] = foam_extend_field[cellI];
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          field.boundaryFieldRef()[patchI][faceI] = foam_extend_field.boundaryField()[patchI][faceI];
        }
      }
    }
    // volVectorField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volVectorField,name)) {
      const Foam_extend_::volVectorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volVectorField,name));
      volVectorField& field(const_cast<volVectorField&>(meshLookupObject(mesh_,volVectorField,name)));
      forAll(field, cellI) {
        forAll(field[cellI], i) {
          field[cellI][i] = foam_extend_field[cellI][i];
        }
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          forAll(field.boundaryField()[patchI][faceI], i) {
            field.boundaryFieldRef()[patchI][faceI][i] = foam_extend_field.boundaryField()[patchI][faceI][i];
          }
        }
      }
    }
    // volTensorField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volTensorField,name)) {
      const Foam_extend_::volTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volTensorField,name));
      volTensorField& field(const_cast<volTensorField&>(meshLookupObject(mesh_,volTensorField,name)));
      forAll(field, cellI) {
        forAll(field[cellI], i) {
          field[cellI][i] = foam_extend_field[cellI][i];
        }
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          forAll(field.boundaryField()[patchI][faceI], i) {
            field.boundaryFieldRef()[patchI][faceI][i] = foam_extend_field.boundaryField()[patchI][faceI][i];
          }
        }
      }
    }
    // volSymmTensorField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::volSymmTensorField,name)) {
      const Foam_extend_::volSymmTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::volSymmTensorField,name));
      volTensorField& field(const_cast<volTensorField&>(meshLookupObject(mesh_,volTensorField,name)));
      forAll(field, cellI) {
        field[cellI].xx()=foam_extend_field[cellI].xx();
        field[cellI].xy()=foam_extend_field[cellI].xy();
        field[cellI].xz()=foam_extend_field[cellI].xz();
        field[cellI].yx()=foam_extend_field[cellI].xy();
        field[cellI].yy()=foam_extend_field[cellI].yy();
        field[cellI].yz()=foam_extend_field[cellI].yz();
        field[cellI].zx()=foam_extend_field[cellI].xz();
        field[cellI].zy()=foam_extend_field[cellI].yz();
        field[cellI].zz()=foam_extend_field[cellI].zz();
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          field.boundaryFieldRef()[patchI][faceI].xx()=foam_extend_field.boundaryField()[patchI][faceI].xx();
          field.boundaryFieldRef()[patchI][faceI].xy()=foam_extend_field.boundaryField()[patchI][faceI].xy();
          field.boundaryFieldRef()[patchI][faceI].xz()=foam_extend_field.boundaryField()[patchI][faceI].xz();
          field.boundaryFieldRef()[patchI][faceI].yx()=foam_extend_field.boundaryField()[patchI][faceI].xy();
          field.boundaryFieldRef()[patchI][faceI].yy()=foam_extend_field.boundaryField()[patchI][faceI].yy();
          field.boundaryFieldRef()[patchI][faceI].yz()=foam_extend_field.boundaryField()[patchI][faceI].yz();
          field.boundaryFieldRef()[patchI][faceI].zx()=foam_extend_field.boundaryField()[patchI][faceI].xz();
          field.boundaryFieldRef()[patchI][faceI].zy()=foam_extend_field.boundaryField()[patchI][faceI].yz();
          field.boundaryFieldRef()[patchI][faceI].zz()=foam_extend_field.boundaryField()[patchI][faceI].zz();
        }
      }
    }
    // surfaceScalarField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceScalarField,name)) {
      const Foam_extend_::surfaceScalarField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceScalarField,name));
      surfaceScalarField& field(const_cast<surfaceScalarField&>(meshLookupObject(mesh_,surfaceScalarField,name)));
      forAll(field, cellI) {
        field[cellI] = foam_extend_field[cellI];
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          field.boundaryFieldRef()[patchI][faceI] = foam_extend_field.boundaryField()[patchI][faceI];
        }
      }
    }
    // surfaceVectorField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceVectorField,name)) {
      const Foam_extend_::surfaceVectorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceVectorField,name));
      surfaceVectorField& field(const_cast<surfaceVectorField&>(meshLookupObject(mesh_,surfaceVectorField,name)));
      forAll(field, cellI) {
        forAll(field[cellI], i) {
          field[cellI][i] = foam_extend_field[cellI][i];
        }
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          forAll(field.boundaryField()[patchI][faceI], i) {
            field.boundaryFieldRef()[patchI][faceI][i] = foam_extend_field.boundaryField()[patchI][faceI][i];
          }
        }
      }
    }
    // surfaceTensorField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceTensorField,name)) {
      const Foam_extend_::surfaceTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceTensorField,name));
      surfaceTensorField& field(const_cast<surfaceTensorField&>(meshLookupObject(mesh_,surfaceTensorField,name)));
      forAll(field, cellI) {
        forAll(field[cellI], i) {
          field[cellI][i] = foam_extend_field[cellI][i];
        }
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          forAll(field.boundaryField()[patchI][faceI], i) {
            field.boundaryFieldRef()[patchI][faceI][i] = foam_extend_field.boundaryField()[patchI][faceI][i];
          }
        }
      }
    }
    // surfaceSymmTensorField
    if (meshFoundObject(foam_extend_mesh_,Foam_extend_::surfaceSymmTensorField,name)) {
      const Foam_extend_::surfaceSymmTensorField& foam_extend_field(meshLookupObject(foam_extend_mesh_,Foam_extend_::surfaceSymmTensorField,name));
      surfaceTensorField& field(const_cast<surfaceTensorField&>(meshLookupObject(mesh_,surfaceTensorField,name)));
      forAll(field, cellI) {
        forAll(field[cellI], i) {
          field[cellI][i] = foam_extend_field[cellI][i];
        }
      }
      forAll(field.boundaryField(), patchI) {
        forAll(field.boundaryField()[patchI], faceI) {
          forAll(field.boundaryField()[patchI][faceI], i) {
            field.boundaryFieldRef()[patchI][faceI][i] = foam_extend_field.boundaryField()[patchI][faceI][i];
          }
        }
      }
    }
  }
}
#endif

// ************************************************************************* //

