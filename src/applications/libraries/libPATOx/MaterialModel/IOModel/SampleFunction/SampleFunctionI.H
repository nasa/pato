/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::SampleFunction::updateProbing()
{
  if(dynamicMesh_) {
    if (firstIteration_) {
      mesh_.oldCellCentres();
    }
    meshSearchPtr_.reset(new meshSearch(mesh_));
    functionPtr_.reset(new functionObjectListTest(runTime, functionDict_));
    functionPtr_->start();
  }

  functionPtr_().executeScalingAll
  (
      scalarMasterFields_ref,
      vectorMasterFields_ref,
      sphericalTensorMasterFields_ref,
      symmTensorMasterFields_ref,
      tensorMasterFields_ref
  ) ;

  // Verify the probes are still in the mesh after moving
  if(dynamicMesh_) {
    forAll(subProbingDictTOC_, dictI) {

      bool patchSampling = false;
      if (subProbingDictTOC_[dictI] == "surfacePatch") {
        patchSampling = true;
      }
      if(!patchSampling) {
        forAll(probingPoints_[dictI], probeI) {
          if (!point_in_mesh(probingPoints_[dictI][probeI])) {
            forAll(probingFields_scalar_,fieldI) {
              if (scalarMasterFields_ref()[fieldI][dictI].size()>0) {
                Field<scalar>& scalarMasterFields_ = scalarMasterFields_ref()[fieldI][dictI];

                Field<scalar> scalarMasterFields_copy(scalarMasterFields_);
                scalarMasterFields_.resize(scalarMasterFields_copy.size()+1);
                forAll(scalarMasterFields_, i) {
                  if(i<probeI) {
                    scalarMasterFields_[i] = scalarMasterFields_copy[i];
                  }
                  if (i == probeI) {
                    scalarMasterFields_[i] = 0;
                  }
                  if(i>probeI) {
                    scalarMasterFields_[i] = scalarMasterFields_copy[i-1];
                  }
                }
              }
            }

            forAll(probingFields_vector_,fieldI) {
              if (vectorMasterFields_ref()[fieldI][dictI].size()>0) {

                Field<vector>& vectorMasterFields_ = vectorMasterFields_ref()[fieldI][dictI];
                Field<vector> vectorMasterFields_copy(vectorMasterFields_);
                vectorMasterFields_.resize(vectorMasterFields_copy.size()+1);
                forAll(vectorMasterFields_, i) {
                  if(i<probeI) {
                    vectorMasterFields_[i] = vectorMasterFields_copy[i];
                  }
                  if (i == probeI) {
                    vectorMasterFields_[i] = vector(0,0,0);
                  }
                  if(i>probeI) {
                    vectorMasterFields_[i] = vectorMasterFields_copy[i-1];
                  }
                }
              }
            }

            forAll(probingFields_tensor_,fieldI) {
              if (tensorMasterFields_ref()[fieldI][dictI].size() > 0) {
                Field<tensor>& tensorMasterFields_ = tensorMasterFields_ref()[fieldI][dictI];
                Field<tensor> tensorMasterFields_copy(tensorMasterFields_);
                tensorMasterFields_.resize(tensorMasterFields_copy.size()+1);
                forAll(tensorMasterFields_, i) {
                  if(i<probeI) {
                    tensorMasterFields_[i] = tensorMasterFields_copy[i];
                  }
                  if (i == probeI) {
                    tensorMasterFields_[i] = tensor(0,0,0,0,0,0,0,0,0);
                  }
                  if(i>probeI) {
                    tensorMasterFields_[i] = tensorMasterFields_copy[i-1];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

inline void Foam::SampleFunction::writeOutput()
{
  updateProbing();

  // Write in the output folders
  if (Pstream::master()) {
    forAll(subProbingDictTOC_, dictI) {
      forAll(probingFields_scalar_, fieldI) {
        os_out_scalar_[dictI][fieldI] << mesh_.time().value() << " ";
        forAll(probingPoints_[dictI], pointI) {
          scalar sample_ = scalarMasterFields_ref()[fieldI][dictI][pointI] ;
          os_out_scalar_[dictI][fieldI] << sample_ << " ";
        }
        os_out_scalar_[dictI][fieldI] << endl;
      }
      forAll(probingFields_vector_, fieldI) {
        for(int coordI =0; coordI < 3; coordI++) {
          os_out_vector_[dictI][fieldI][coordI] << mesh_.time().value() << " ";
          forAll(probingPoints_[dictI], pointI) {
            scalar sample_ = vectorMasterFields_ref()[fieldI][dictI][pointI][coordI];
            os_out_vector_[dictI][fieldI][coordI] << sample_ << " ";
          }
          os_out_vector_[dictI][fieldI][coordI] << endl;
        }
      }
      forAll(probingFields_tensor_, fieldI) {
        for(int coordI =0; coordI < 9; coordI++) {
          os_out_tensor_[dictI][fieldI][coordI] << mesh_.time().value() << " ";
          forAll(probingPoints_[dictI], pointI) {
            scalar sample_ = tensorMasterFields_ref()[fieldI][dictI][pointI][coordI];
            os_out_tensor_[dictI][fieldI][coordI] << sample_ << " ";
          }
          os_out_tensor_[dictI][fieldI][coordI] << endl;
        }
      }
    }
  }
}


//- list of the points to probe per sub dictionary in the probing dictionary
inline Foam::PtrList<Foam::List<Foam::vector> >& Foam::SampleFunction::probingPoints()
{
  return probingPoints_;
}

//- list of the scalar fields to probe
inline Foam::PtrList<Foam::word>& Foam::SampleFunction::probingFields_scalar()
{
  return probingFields_scalar_;
}

//- list of the vector fields to probe
inline Foam::PtrList<Foam::word>& Foam::SampleFunction::probingFields_vector()
{
  return probingFields_vector_;
}

//- list of the tensor fields to probe
inline Foam::PtrList<Foam::word>& Foam::SampleFunction::probingFields_tensor()
{
  return probingFields_tensor_;
}

inline Foam::fileName& Foam::SampleFunction::dictPath()
{
  return dictPath_;
}

inline Foam::autoPtr<Foam::List<Foam::List<Foam::Field<Foam::scalar> > > >& Foam::SampleFunction::scalarMasterFields()
{
  return scalarMasterFields_ref;
}

inline void Foam::SampleFunction::lower_string(string& str)
{
  for(int i=0; str[i]!='\0'; i++) {
    if (str[i] >= 'A' && str[i] <= 'Z')    //checking for uppercase characters
      str[i] = str[i] + 32;         //converting uppercase to lowercase
  }
}

inline void Foam::SampleFunction::changeUpperLowercaseName(PtrList<word>& list)
{
  List<word> lower_list(list.size());
  forAll(lower_list, i) {
    lower_list[i]=list[i];
    lower_string(lower_list[i]);
  }
  List<int> count(list.size(),0);
  forAll(lower_list, i) {
    forAll(lower_list, j) {
      if (lower_list[i] == lower_list[j] && i!=j) {
        count[j]++;
        lower_list[j]+="_";
      }
    }
  }
  forAll(list, i) {
    for (int j=0; j<count[i]; j++) {
      list[i]+="_";
    }
  }
}

inline bool Foam::SampleFunction::point_in_mesh(const vector& point)
{
  if (Pstream::parRun()) {
    List<int> list_find_cell(Pstream::nProcs(),0);
    list_find_cell[Pstream::myProcNo()]=meshSearchPtr_().findCell(point);
    reduce(list_find_cell, sumOp<List<int> >());
    forAll(list_find_cell, i) {
      if (list_find_cell[i]>=0) {
        return true;
      }
    }
    return false;
  } else {
    return (meshSearchPtr_().findCell(point)>=0);
  }
}

// ************************************************************************* //
