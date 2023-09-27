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

#include "TecplotReader.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TecplotReader::TecplotReader
(
    const fileName& file
)
  :
fileName_(changeEnviVar(file)),
zoneI(0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TecplotReader::~TecplotReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::TecplotReader::readHeader(IFstream& dataFile_)
{
  int index_=0; // index 0 = title, index 1 = variables

  while(dataFile_) {
    string line;
    dataFile_.getLine(line);
    if (line.empty()||onlySpaces(line)) { // empty line or only space/tab
      break;
    }
    if (index_==0) { // TITLE
      title_ = line;
      index_=1;
      continue;
    }
    if (index_ == 1) { // VARIABLES
      int zone_t_found  = line.find("ZONE T");
      int variables_found  = line.find("VARIABLES");
      if(zone_t_found<0) { // ZONE T not found
        if(variables_found<0) { // VARIABLES not found
          variables_.append(line.replaceAll("\"","")); // add variables
        } else {  // to get the first variable
          line.replaceAll("VARIABLES","");
          line.replaceAll("=","");
          line.replaceAll(" ","");
          line.replaceAll("\"","");
          if (!(line.empty()||onlySpaces(line))) {
            variables_.append(line);
          }
        }
      } else {
        data_.resize(variables_.size()); // resize the data in function of the variables
        break;
      }
    }
  }
}

void Foam::TecplotReader::readZoneHeader(IFstream& dataFile_)
{
  int index_ = 0;

  while(dataFile_) {
    string line;
    dataFile_.getLine(line);
    if (line.empty()||onlySpaces(line)) { // empty line or only space/tab
      break;
    }
    if (index_==0) { // STRANDID || SOLUTIONTIME
    }
    if(index_==1) { // I,J,K or Elements, Nodes

      int foundFETriangle= line.find("FETriangle");
      int foundOrdered= line.find("Ordered");

      if (foundFETriangle>=0) {
        zoneType.append("FETriangle");
      }
      if (foundOrdered>=0) {
        zoneType.append("Ordered");
      }

      if (foundFETriangle<0 && foundOrdered<0) {
        FatalErrorInFunction << "readTecplot not implemented for " << line << exit(FatalError);
      }

      if (foundOrdered>=0) { // I,J,K
        I_zones_.append(findScalarBetweenStringsAndStrips(line,"I=", ","));
        J_zones_.append(findScalarBetweenStringsAndStrips(line,"J=", ","));
        K_zones_.append(findScalarBetweenStringsAndStrips(line,"K=", ","));
      }
    }
    if (index_==2) { // DATAPACKING
      int loc_ = line.find("BLOCK");
      if (loc_ < 0) {
        FatalErrorInFunction << "Only \"DATAPACKING=BLOCK\" is implemented" << exit(FatalError);

      }
    }
    if (index_==3) { // DT
      break;
    }
    index_++;
  }

}


void Foam::TecplotReader::readZoneData(IFstream& dataFile_)
{
//    ijk_nodes_.resize(zoneI+1);
//    xyz_nodes_.resize(zoneI+1);
  int current_variable_ = 0;
  int current_points_number_=0;
  int max_points_ = I_zones_[zoneI]*J_zones_[zoneI]*K_zones_[zoneI];
//  int i_ = 0;
//  int j_ = 0;
//  int k_ = 0;

  while(dataFile_) { // ZONES DATA: VAR_1_I VAR_1_J VAR_1_K VAR_2_I ...
    string line;
    dataFile_.getLine(line);
    if (line.empty()||onlySpaces(line)) { // empty line or only space/tab
      break;
    }
    int zone_t_found  = line.find("ZONE T");
    if (zone_t_found>=0) { // STARTS over if ZONE T is found
      break;
    }
    lineStream_.reset(new IStringStream(line));
    if(zoneType[zoneI]=="Ordered") { // Ordered type

      // DATAPACKING = BLOCK:
      // Read data: x_0 ... x_max_points_ y_0 ... y_max_points ...
      scalar value_;
      while(lineStream_()) {
        lineStream_() >> value_;
        if (current_variable_>=data_.size()) {
          FatalErrorInFunction << "Problem reading the Tecplot file" << exit(FatalError);
        }
        data_[current_variable_].append(value_);
//        if (current_variable_==0){
//            vector ijk(i_,j_,k_);
//            ijk_nodes_[zoneI].append(ijk);
//            i_++;
//            if (i_>I_zones_[zoneI]-1){
//                j_++;
//                i_=0;
//            }
//            if (j_ > J_zones_[zoneI]-1){
//                k_++;
//                j_=0;
//                i_=0;
//            }
//        }
//        if(current_variable_<3){
//            xyz_nodes_[zoneI][current_variable_].append(value_);
//        }
        current_points_number_++;
        if(current_points_number_>= max_points_) {
          current_variable_++;
          current_points_number_=0;
        }
      }
    }

  }
}

void Foam::TecplotReader::writePlot3D(word name, List<scalarList>& data)
{
  int nZones = I_zones_.size();
  if ((J_zones_.size()!=nZones)||(K_zones_.size()!=nZones)) {
    FatalErrorInFunction << "I,J,K does not have the same number of zones." << exit(FatalError);
  }
  OFstream os_out(name+"_ASCII.g");
  Info << "Writing " << name +"_ASCII.g" << endl;

  OFstream * os_out2 = new OFstream(name+"_ASCII.f");;
  Switch writeF=false;
  if (data.size()>3) {
    Info << "Writing " << name +"_ASCII.f" << endl;
    writeF=true;
  }

  os_out << nZones << endl;

  forAll(I_zones_, zI) {
    os_out << I_zones_[zI] << " " << J_zones_[zI] << " " << K_zones_[zI] << endl;
  }
  List<int> index(3);
  forAll(index, i) {
    index[i]=0;
  }
  forAll(I_zones_, zI) {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < I_zones_[zI]*J_zones_[zI]*K_zones_[zI]; j++) {
        os_out << data[i][index[i]] << " ";
        index[i]++;
      }
    }
    os_out << endl;
  }

  if(writeF) {
    int sizeFunctions = data.size()-3;
    *os_out2 << nZones << endl;
    forAll(I_zones_, zI) {
      *os_out2  << I_zones_[zI] << " " << J_zones_[zI] << " " << K_zones_[zI] << " " << sizeFunctions << endl;
    }
    List<int> index(sizeFunctions);
    forAll(index, i) {
      index[i]=0;
    }
    forAll(I_zones_, zI) {
      for(int i = 0; i < sizeFunctions; i++) {
        for(int j = 0; j < I_zones_[zI]*J_zones_[zI]*K_zones_[zI]; j++) {
          *os_out2 << data[i+3][index[i]] << " ";
          index[i]++;
        }
      }
      *os_out2 << endl;
    }
  }
}

void Foam::TecplotReader::splitData(const List<scalarList>& data)
{
  ijk_nodes_.resize(I_zones_.size());
  xyz_nodes_.resize(I_zones_.size());
  zoneI=0;
  int currentPoint=0;
  int i_=0;
  int j_=0;
  int k_=0;
  forAll(data[0], nI) {
    int maxPoints=I_zones_[zoneI]*J_zones_[zoneI]*K_zones_[zoneI];
    if (currentPoint>=maxPoints) {
      zoneI++;
      currentPoint=0;
      i_=0;
      j_=0;
      k_=0;
    }
    vector ijk(i_,j_,k_);
    ijk_nodes_[zoneI].append(ijk);
    i_++;
    if (i_>I_zones_[zoneI]-1) {
      j_++;
      i_=0;
    }
    if (j_ > J_zones_[zoneI]-1) {
      k_++;
      j_=0;
      i_=0;
    }

    xyz_nodes_[zoneI].append(vector(data[0][nI],data[1][nI],data[2][nI]));
    currentPoint++;
  }

}

void Foam::TecplotReader::readFile()
{


  IFstream dataFile_(fileName_);
  if (!dataFile_.good()) {
    FatalErrorInFunction
        << "Cannot read file " << fileName_
        << exit(FatalError);
  }

  readHeader(dataFile_);
  while(dataFile_) {
    readZoneHeader(dataFile_);
    readZoneData(dataFile_);
    zoneI++;
  }
  nZones_=zoneI;
}

int Foam::TecplotReader::nZones()
{
  return nZones_;
}

scalarList Foam::TecplotReader::I_zones()
{
  return I_zones_;
}

scalarList Foam::TecplotReader::J_zones()
{
  return J_zones_;
}

scalarList Foam::TecplotReader::K_zones()
{
  return K_zones_;
}

const List<scalarList>& Foam::TecplotReader::data()
{
  return data_;
}

const List<List<vector>>& Foam::TecplotReader::ijk_nodes()
{
  return ijk_nodes_;
}

const List<List<vector>>& Foam::TecplotReader::xyz_nodes()
{
  return xyz_nodes_;
}

// ************************************************************************* //
