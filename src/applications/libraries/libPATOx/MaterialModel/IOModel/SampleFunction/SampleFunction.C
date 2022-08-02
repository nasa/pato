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

#include "SampleFunction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SampleFunction::SampleFunction
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName,
    const word& probingDictName
):
mesh_(mesh),
runTime(mesh.time()),
probingDictName_(probingDictName),
dictPath_(mesh.time().system() + "/" + dictName + "/"+ probingDictName_),
C_initial(mesh.C()),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
functionDict_(dictionary::null),
rank_(Pstream::myProcNo()),
firstIteration_(true)
{
  IOdictionary probingDict_
  (
      IOobject
      (
          probingDictName_,
          mesh.time().system(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  //- list of the fields to probe
  wordList tmp_list(probingDict_.lookup("fields"));

  // Verify the fields are stored in mesh
  foundFieldsInMeshAll(mesh, tmp_list); // includes "all" field

  fileName outputDirectory_ = "output/"+dictName;
  forAll(tmp_list, fieldI) {
    if (mesh_.objectRegistry::foundObject<volScalarField>(tmp_list[fieldI])) {
      if (!isDir(outputDirectory_+"/scalar")) {
        system("mkdir -p " + outputDirectory_ + "/scalar");
      }
      probingFields_scalar_.append(new word(tmp_list[fieldI]));
    }
    if (mesh_.objectRegistry::foundObject<volVectorField>(tmp_list[fieldI])) {
      if (!isDir(outputDirectory_+"/vector")) {
        system("mkdir -p " + outputDirectory_ + "/vector");
      }
      probingFields_vector_.append(new word(tmp_list[fieldI]));
    }
    if (mesh_.objectRegistry::foundObject<volTensorField>(tmp_list[fieldI])) {
      if (!isDir(outputDirectory_+"/tensor")) {
        system("mkdir -p " + outputDirectory_ + "/tensor");
      }
      probingFields_tensor_.append(new word(tmp_list[fieldI]));
    }
  }

  // Change the name of the element list if the name already exists in lower/uppercase
  changeUpperLowercaseName(probingFields_scalar_);
  changeUpperLowercaseName(probingFields_vector_);
  changeUpperLowercaseName(probingFields_tensor_);

  const word typeDict_(probingDict_.lookup("type"));

  dictionary subProbingDict_;
  ITstream itstream_(probingDict_.lookup(typeDict_));
  token token_;
  itstream_ >> token_ >> subProbingDict_;

  wordList subProbingDictTOC_tmp(subProbingDict_.toc());
  subProbingDictTOC_.resize(subProbingDictTOC_tmp.size());
  forAll(subProbingDictTOC_tmp, wordI) {
    subProbingDictTOC_.set(wordI, new word(subProbingDictTOC_tmp[wordI]));
  }

  forAll(subProbingDictTOC_, wordI) {
    listSubProbingDict_.append(subProbingDict_.subDict(subProbingDictTOC_[wordI]));
  }

  probingPoints_.resize(subProbingDictTOC_.size());
  forAll(subProbingDictTOC_, dictI) {
    probingPoints_.set
    (
        dictI,
        new List<vector>
        (
            listSubProbingDict_[dictI].lookup("points")
        )
    );
  }

  forAll(subProbingDictTOC_, dictI) {
    os_out_scalar_.append(new PtrList<OFstream>);
    forAll(probingFields_scalar_, fieldI) {
      os_out_scalar_[dictI].append(new OFstream(outputDirectory_+"/scalar/"+probingFields_scalar_[fieldI]+"_"+subProbingDictTOC_[dictI]));
      os_out_scalar_[dictI][fieldI] << "//t(s) ";
      forAll(probingPoints_[dictI], pointI) {
        os_out_scalar_[dictI][fieldI] << "probe" << pointI  << "("<< probingPoints_[dictI][pointI][0] << ","
                                      << probingPoints_[dictI][pointI][1] << "," << probingPoints_[dictI][pointI][2] << ") ";
      }
      os_out_scalar_[dictI][fieldI] << endl;
    }
  }

  wordList nameVector;
  nameVector.append("x");
  nameVector.append("y");
  nameVector.append("z");

  forAll(subProbingDictTOC_, dictI) {
    os_out_vector_.append(new PtrList<PtrList<OFstream> >);
    forAll(probingFields_vector_, fieldI) {
      os_out_vector_[dictI].append(new PtrList<OFstream>);
      forAll(nameVector, coordI) {
        os_out_vector_[dictI][fieldI].append(new OFstream(outputDirectory_+"/vector/"+probingFields_vector_[fieldI]+nameVector[coordI]+"_"+subProbingDictTOC_[dictI]));
        os_out_vector_[dictI][fieldI][coordI] << "//t(s) ";
        forAll(probingPoints_[dictI], pointI) {
          os_out_vector_[dictI][fieldI][coordI] << "probe" << pointI  << "("<< probingPoints_[dictI][pointI][0] << ","
                                                << probingPoints_[dictI][pointI][1] << "," << probingPoints_[dictI][pointI][2] << ") ";
        }
        os_out_vector_[dictI][fieldI][coordI] << endl;
      }
    }
  }

  wordList nameTensor;
  nameTensor.append("xx");
  nameTensor.append("xy");
  nameTensor.append("xz");
  nameTensor.append("yx");
  nameTensor.append("yy");
  nameTensor.append("yz");
  nameTensor.append("zx");
  nameTensor.append("zy");
  nameTensor.append("zz");

  forAll(subProbingDictTOC_, dictI) {
    os_out_tensor_.append(new PtrList<PtrList<OFstream> >);
    forAll(probingFields_tensor_, fieldI) {
      os_out_tensor_[dictI].append(new PtrList<OFstream>);
      forAll(nameTensor, coordI) {
        os_out_tensor_[dictI][fieldI].append(new OFstream(outputDirectory_+"/tensor/"+probingFields_tensor_[fieldI]+nameTensor[coordI]+"_"+subProbingDictTOC_[dictI]));
        os_out_tensor_[dictI][fieldI][coordI] << "//t(s) ";
        forAll(probingPoints_[dictI], pointI) {
          os_out_tensor_[dictI][fieldI][coordI] << "probe" << pointI  << "("<< probingPoints_[dictI][pointI][0] << ","
                                                << probingPoints_[dictI][pointI][1] << "," << probingPoints_[dictI][pointI][2] << ") ";
        }
        os_out_tensor_[dictI][fieldI][coordI] << endl;
      }
    }
  }

  if (dictName != word::null) {
    probingDict_.set("region", dictName);
  }
  dictionary functionDict1_ = dictionary::null;
  functionDict1_.add
  (
      dictionaryEntry(probingDictName_, functionDict1_, probingDict_)
  );

  functionDict_.add
  (
      dictionaryEntry("functions", functionDict_, functionDict1_)
  );

  if(!dynamicMesh_) {
    meshSearchPtr_.reset(new meshSearch(mesh_));
    functionPtr_.reset(new functionObjectListTest(runTime, functionDict_));
    functionPtr_->start();
  }

  // Write first time step
  writeOutput();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SampleFunction::~SampleFunction()
{}


// ************************************************************************* //
