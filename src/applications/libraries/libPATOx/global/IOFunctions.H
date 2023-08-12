/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    IO functions

Details

\*---------------------------------------------------------------------------*/

#ifndef IOFunctionsFoam_H
#define IOFunctionsFoam_H

#include "fvCFD.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Return the global face center vectors
inline List<vector> globalFaceCenters(const fvMesh& mesh, const word& patchName, const word& regionName)
{
  const label patchID = mesh.boundaryMesh().findPatchID(patchName); // patch id
  if(patchID<0) {
    FatalErrorInFunction << "patch " << patchName << " not found in mesh." << exit(FatalError);
  }

  // If single processor, return mesh.Cf().bF()[patchID]
  if (!Pstream::parRun()) {
    return mesh.Cf().boundaryField()[patchID];
  }

  const int nProc = Pstream::nProcs(); // number of processors
  const int rank = Pstream::myProcNo(); // rank of processors

  // Parallel boundary info
  polyBoundaryMesh boundaryIO
  (
      IOobject
      (
          "boundary",
          "constant/"+regionName+"/polyMesh",
          mesh.time(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      ),
      mesh
  );
  wordList boundary_names=boundaryIO.names();
  label boundary_name_id=-1;
  forAll(boundary_names, i) {
    if (boundary_names[i]==patchName) {
      boundary_name_id = i;
      break;
    }
  }
  if (boundary_name_id<0) {
    FatalError << patchName << " not found in mesh boundary." << exit(FatalError);
  }
  int nFaces = boundaryIO[boundary_name_id].size();
  label startFace = boundaryIO[boundary_name_id].start();

  // Local faces processor address
  word local_face_path=word("constant/"+regionName+"/polyMesh/faceProcAddressing");
  labelIOList localFaceProcAddr
  (
      IOobject
      (
          local_face_path,
          mesh.time(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );
  List<label> localFacePatch(nFaces);
  forAll(localFacePatch, i) {
    localFacePatch[i]=localFaceProcAddr[startFace+i];
  }

  // Compute the maximum face id
  label localFacePatch_max=0;
  if (localFacePatch.size()>0) {
    localFacePatch_max=max(localFacePatch);
  }
  List<label> globalFacePatch_max(nProc);
  forAll(globalFacePatch_max, i) {
    globalFacePatch_max[i]=0;
  }
  globalFacePatch_max[rank]=localFacePatch_max;
  reduce(globalFacePatch_max, sumOp<List<label> >());

  // Compute the minimum face id
  List<label> localFacePatch_pos;
  forAll(localFacePatch, i) {
    if(localFacePatch[i]>=0) {
      localFacePatch_pos.append(localFacePatch[i]);
    }
  }
  label localFacePatch_min=0;
  if (localFacePatch_pos.size()>0) {
    localFacePatch_min=min(localFacePatch_pos);
  }
  List<label> globalFacePatch_min(nProc);
  forAll(globalFacePatch_min, i) {
    globalFacePatch_min[i]=0;
  }
  globalFacePatch_min[rank]=localFacePatch_min;
  reduce(globalFacePatch_min, sumOp<List<label> >());
  List<label> globalFacePatch_min_non_zero;
  forAll(globalFacePatch_min, i) {
    if (globalFacePatch_min[i]>0) {
      globalFacePatch_min_non_zero.append(globalFacePatch_min[i]);
    }
  }

  // Global face centers of patch "patchName"
  int max_globalFacePatch = max(globalFacePatch_max);
  int min_globalFacePatch = min(globalFacePatch_min_non_zero);
  int global_size_Cf = max_globalFacePatch - min_globalFacePatch + 1;
  List<vector> global_Cf(global_size_Cf); // global face centers of patch "patchName"
  forAll(global_Cf, i) {
    global_Cf[i]=vector(0,0,0);
  }
  forAll(localFacePatch, i) {
    if (localFacePatch[i]>=0) {
      global_Cf[localFacePatch[i]-min_globalFacePatch]=mesh.Cf().boundaryField()[patchID][i];
    }
  }
  reduce(global_Cf, sumOp<List<vector> >());
  return global_Cf;
}


//- remove C++ comments from a list of string
inline List<string> removeComments(List<string> lines)
{
  bool comment_start=false;
  int found_start;
  int found_end;

  forAll(lines, i) {
    string& line=lines[i];
    line=line.substr(0,line.find("//"));
    found_start = line.find("/*");
    found_end = line.find("*/");
    if (found_start >= 0) {
      comment_start=true;
    }
    if (comment_start) {
      if (found_start >= 0) {
        if (found_end >= 0) {
          if (found_end>found_start) {
            line=line.erase(found_start, found_end-found_start+2);
            comment_start=false;
          }
        } else {
          line=line.erase(found_start, line.size()-found_start);
        }
      } else {
        if (found_end >= 0) {
          line=line.erase(0, found_end+2);
          comment_start=false;
        } else {
          line="";
        }
      }
    }
  }
  return lines;
}

//- return true if the string has only tab or spaces
inline bool onlySpaces(string testString)
{
  string newString = testString;
  newString.replaceAll(" ", "");
  newString.replaceAll("\t", "");
  return newString=="";
}

//- Trim a string (remove the space before and after)
inline string trim_space(const string& str)
{
  size_t first = str.find_first_not_of(' ');
  if (string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

//- Trim a string (remove the tab before and after)
inline string trim_tab(const string& str)
{
  size_t first = str.find_first_not_of('\t');
  if (string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of('\t');
  return str.substr(first, (last - first + 1));
}

//- Trim a string (remove the space and tab before and after)
inline string trim(const string& str)
{
  string new_str=str;
  string previous_str="";
  while(new_str!=previous_str) {
    previous_str=new_str;
    new_str=trim_space(new_str);
    new_str=trim_tab(new_str);
  }
  return new_str;
}

//- Read a new line from a text file
class lineReader
{
 private:
  //- New line from a text file
  string line_;
  //- Words from the split of the new line
  List<string> splitLine_;
  //- Word index
  int wordIndex_;
 public:
  inline void readNewLine(string line);
  inline void readNext(double& value);
  inline void readNext(string& value);
  inline List<string> splitStringBySpace(string str);
  inline string line() {
    return line_;
  };
  inline List<string> splitLine() {
    return splitLine_;
  };
  inline int wordIndex() {
    return wordIndex_;
  };
};

//- Read a new line
inline void lineReader::readNewLine(string line)
{
  line_=line;
  splitLine_=splitStringBySpace(line);
  wordIndex_=0;
}

//- Read the next numerical value from the new line
inline void lineReader::readNext(double& value)
{
  if (wordIndex_<splitLine_.size()) {
    value = stof(splitLine_[wordIndex_]);
    wordIndex_++;
  }
}

//- Read the next string value from the new line
inline void lineReader::readNext(string& value)
{
  if (wordIndex_<splitLine_.size()) {
    value = splitLine_[wordIndex_];
    wordIndex_++;
  }
}

//- Split a string by space
inline List<string> lineReader::splitStringBySpace(string str)
{
  List<string> words;
  string word = "";
  for (unsigned int i = 0; i < str.size(); i++) {
    char x=str[i];
    if (x == ' ') {
      if (word!="") {
        words.append(word);
      }
      word = "";
    } else {
      word = word + x;
    }
  }
  if (word!="") {
    words.append(word);
  }
  return words;
}

//- Returns true if elem is found in list
template <typename Type> inline  bool foundInList(Type elem, List<Type> list)
{
  bool found = false;
  forAll(list, elemI) {
    if(list[elemI]==elem) {
      found = true;
      break;
    }
  }
  return found;
}

//- Returns the index of the first elem in list
template <typename Type> inline int indexInList(Type elem, List<Type> list)
{
  if (!foundInList<Type>(elem,list)) {
    FatalError << elem << " not found in " << list << exit(FatalError);
  }
  int index = -1;
  forAll(list, elemI) {
    if(list[elemI]==elem) {
      index = elemI;
      break;
    }
  }
  return index;
}

// Splits a name to catch exponential
inline Foam::wordList splitName
(
    const word& name,
    const int nCmpt
)
{
  wordList cmpts(nCmpt);

  string::size_type beg=0, end=0, endb=0, endc=0;
  int i = 0;

  while
  (
      (endb = name.find('e', beg)) != string::npos
      || (endc = name.find('E', beg)) != string::npos
  ) {
    if (endb == string::npos) {
      end = endc;
    } else if ((endc = name.find('E', beg)) != string::npos) {
      end = std::min(endb, endc);
    } else {
      end = endb;
    }

    if (beg < end) {
      cmpts[i] = name.substr(beg, end-beg);
      i++;

      // If the number of number of components in the name
      // is greater than nCmpt return an empty list
      if (i == nCmpt) {
        wordList list;
        return list;
      }
    }
    beg = end + 1;
  }

  // If the number of number of components in the name is not equal to nCmpt
  // return an empty list
  if (i + 1 != nCmpt) {
    wordList list;
    return list;
  }

  if (beg < name.size()) {
    cmpts[i] = name.substr(beg, string::npos);
  }

  return cmpts;
}

// Return true if string is a number
inline bool isNumber(const std::string& s)
{
  try {
    wordList splitNames = splitName(s,2);
    scalar numAfterExp = 0;

    if (splitNames.size()>1) {
      numAfterExp = stof(splitName(s,2)[1]);
    }
    if ( numAfterExp<37 && numAfterExp>-37) { // float limit for the exponential
      stof(s);
    } else {
      return false;
    }
  } catch(...) {
    return false;
  }
  return true;
}

//- Return the number between two strings
inline scalar findScalarBetweenStrings(const string& line, const string& str_start, const string& str_end)
{
  int first = line.find_first_of(str_start);
  int last = line.find_first_of(str_end);
  if (first < 0 || last < 0) {
    FatalErrorInFunction << str_start << " or " << str_end << " not found in " << line << exit(FatalError);
  }
  string num_ = line.substr(first+str_start.size(),last-(first+str_start.size()));
  if (isNumber(num_)) {
    scalar number_=stof(num_);
    return number_;
  } else {
    FatalErrorInFunction << "No scalar found between " << str_start << " and " << str_end << " (first found) in " << line << exit(FatalError);
  }

  return 0;
}

//- Return the number between two strings and strips the line
inline scalar findScalarBetweenStringsAndStrips(string& line, const string& str_start, const string& str_end)
{
  int first = line.find_first_of(str_start);
  int last = line.find_first_of(str_end);
  if (first < 0 || last < 0) {
    FatalErrorInFunction << str_start << " or " << str_end << " not found in " << line << exit(FatalError);
  }
  string num_ = line.substr(first+str_start.size(),last-(first+str_start.size()));
  if (isNumber(num_)) {
    scalar number_=stof(num_);
    int strip_= last+str_end.size();
    int size_ = line.size();
    if (strip_ >= 0 && strip_< size_) {
      line = line.substr(strip_,size_);
    }
    return number_;
  } else {
    FatalErrorInFunction << "No scalar found between " << str_start << " and " << str_end << " (first found) in " << line << exit(FatalError);
  }

  return 0;
}

//- Change environment variable in fileName to full path
inline fileName changeEnviVar(const fileName filenameOrigin_)
{
  // Change environment variable from filename_
  wordList enviVariables_;
  fileName fileNameNew_ = filenameOrigin_;
  fileName fileName_ = filenameOrigin_;
  int i = 0;
  while(i < 100) {
    i++;
    int first = fileNameNew_.find("$");
    int last = fileNameNew_.find("/");
    if (first < 0 ) {
      break;
    }
    string strNew = fileNameNew_.substr (first+1,last-first-1);
    enviVariables_.append(strNew);
    fileNameNew_=fileNameNew_.replaceAll( fileNameNew_.substr (first,last-first),"");

  }

  forAll(enviVariables_, enviI) {
    fileName envi_dir = getEnv(enviVariables_[enviI]);
    fileName_.replaceAll("$"+enviVariables_[enviI] , envi_dir);
  }
  return fileName_;
}

//- Read the scalar values in a file and return table of the values
inline List<scalarList> readFileData(const fileName filename_)
{
  fileName fileName_ = changeEnviVar(filename_);
  IFstream dataFile(fileName_);
  if (!dataFile.good()) {
    FatalErrorInFunction
        << "Cannot read file " << filename_
        << exit(FatalError);
  }

  List<scalarList> data;
  List<scalarList> data_t;
  IStringStream * lineStream;
  scalar dataI;
  List<string> lines;
  while (dataFile) {
    string line;
    dataFile.getLine(line);
    lines.append(line);
  }
  lines=removeComments(lines);
  forAll(lines, i) {
    string line=lines[i];
    if (line.empty()||onlySpaces(line)) { // empty line or only space/tab
      continue;
    }
    line=trim(line);
    lineStream = new IStringStream(line);
    scalarList data_line;
    while(*lineStream) {
      *lineStream >> dataI;
      data_line.append(dataI);
    }
    data.append(data_line);
  }

  // transpose the data
  if (data.size()>0) {
    int firstLength = data[0].size();
    forAll(data, rowI) {
      if (data[rowI].size() != firstLength) {
        FatalErrorInFunction << "The number of columns in line " <<  rowI << " from " << filename_ << " are different to the first row.";
        FatalErrorInFunction << exit(FatalError);
      }
    }
    data_t.resize(data[0].size());
    forAll(data_t, lineI) {
      data_t[lineI].resize(data.size());
    }
    forAll(data,rowI) {
      forAll(data[rowI], columnI) {
        data_t[columnI][rowI] = data[rowI][columnI];
      }
    }
  }

  return data_t;
}

// Return the centroid of the input triangles
inline List<vector> centroidTriangle(const List<vector>& nodesElements_, const List<vector>& xyzNodes_)
{
// nodesElements_.size() = nElements, xyzNodes.size() = nNodes
  List<vector> coordCentroid_;
  forAll(nodesElements_,elemI) {
    vector node1_= xyzNodes_[nodesElements_[elemI][0]];
    vector node2_ = xyzNodes_[nodesElements_[elemI][1]];
    vector node3_ = xyzNodes_[nodesElements_[elemI][2]];
    vector centroid_= (node1_+node2_+node3_)/3;
    coordCentroid_.append(centroid_);
  }
  return coordCentroid_;
}


//- Read Tecplot data file
inline List<scalarList> readTecplotFileData(const fileName filename_)
{
  fileName fileName_ = changeEnviVar(filename_);
  IFstream dataFile(fileName_);
  if (!dataFile.good()) {
    FatalErrorInFunction
        << "Cannot read file " << filename_
        << exit(FatalError);
  }

  IStringStream * lineStream;

  int index_=0; // index 0 = title, index 1 = variables, index 2 = zones titles, index 3 = zones data
  int current_line_zones_titles_=0;
  int current_points_number_=0;
  int current_variable_=0;
  int passedLine = 0;
  int currentI_=-1; // Ordered type
  int currentJ_=-1; // Ordered type
  int currentK_=-1; // Ordered type
  int currentElements_=-1; // FETriangle type
  int currentNodes_ = -1; // FETriangle type
  int nElements_=0; // FETriangle type
  int nNodes_=0; // FETriangle type
  bool block_flag = false;
  bool point_flag = false;

  scalarList varlocation_;
  word title_;
  wordList variables_;
  List<scalarList> data_;
  List<vector> nodesElements_;

  // Reading data file
  while(dataFile) {
    string line;
    dataFile.getLine(line);
    if (line.empty()||onlySpaces(line)) { // empty line or only space/tab
      continue;
    }
    if (index_==0) { // TITLE
      title_ = line;
      index_=1;
      continue;
    }
    if (index_ == 1) { // VARIABLES
      int zone_t_found  = line.find("ZONE T");
      int variables_found  = line.find("VARIABLES");
      if(zone_t_found<0) {
        if(variables_found<0) {
          variables_.append(line.replaceAll("\"",""));
        } else {
          line.replaceAll("VARIABLES","");
          line.replaceAll("=","");
          line.replaceAll(" ","");
          line.replaceAll("\"","");
          if (!(line.empty()||onlySpaces(line))) {
            variables_.append(line);
          }
        }
      } else {
        index_=2;
        data_.resize(variables_.size());
        continue;
      }
    }
    if (index_ == 2) { // ZONES TITLES
      if (current_line_zones_titles_==0) { // STRANDID || SOLUTIONTIME
        line=trim(line);
        lineStream = new IStringStream(line);
        int foundFETriangle= line.find("FETriangle");
        int foundOrdered= line.find("Ordered");
        if (foundFETriangle>=0 || foundOrdered>=0) {
          current_line_zones_titles_++;
        }
      }
      if(current_line_zones_titles_==1) {
        line=trim(line);
        lineStream = new IStringStream(line);
        int foundFETriangle= line.find("FETriangle");
        int foundOrdered= line.find("Ordered");

        if (foundFETriangle<0 && foundOrdered<0) {
          FatalErrorInFunction << "readTecplot not implemented for " << line << exit(FatalError);
        }
        if (foundFETriangle>=0) { //Elements, Nodes
          currentNodes_ = findScalarBetweenStringsAndStrips(line,"Nodes=", ",");
          currentElements_ =  findScalarBetweenStringsAndStrips(line,"Elements=", ",");
          nNodes_+=currentNodes_;
          nElements_+=currentElements_;
        }
        if (foundOrdered>=0) { // I,J,K
          currentI_ = findScalarBetweenStringsAndStrips(line,"I=", ",");
          currentJ_ = findScalarBetweenStringsAndStrips(line,"J=", ",");
          currentK_ = findScalarBetweenStringsAndStrips(line,"K=", ",");
        }
      }
      if (current_line_zones_titles_==2) { // DATAPACKING
        int loc_ = line.find("BLOCK");
        block_flag=true;
        if (loc_ < 0) {
          loc_ = line.find("POINT");
          point_flag=true;
          block_flag=false;
          if (loc_ < 0) {
            FatalErrorInFunction << "Only \"DATAPACKING=BLOCK\" or \"DATAPACKING=POINT\" are implemented" << exit(FatalError);
          }
        }
      }
      if (current_line_zones_titles_==3) { // DT
        if (currentElements_>=0 && passedLine==0) { // FETriangle type VARLOCATION
          word varloc_ = line;
          int loc_ = varloc_.find("-");
          if (loc_ < 0) {
            FatalErrorInFunction << "VARLOCATION implemented only for [a-b]=CELLCENTERED" << exit(FatalError);
          }
          int start_ = findScalarBetweenStrings(varloc_,"VARLOCATION=([", "-");
          int end_ = findScalarBetweenStrings(varloc_,"-", "]=CELLCENTERED)");
          for (int i = start_; i<=end_; i++) {
            varlocation_.append(i-1);
          }
          passedLine++;
          continue;
        }
        passedLine=0;
        current_line_zones_titles_=0;
        index_=3;
        continue;
      }
      current_line_zones_titles_++;
    }

    if (index_ == 3) { // ZONES DATA: VAR_1_I VAR_1_J VAR_1_K VAR_2_I ...
      int zone_t_found  = line.find("ZONE T");

      if (zone_t_found>=0) { // STARTS over if ZONE T is found
        current_variable_ = 0;
        current_points_number_=0;
        index_ = 2;
        continue;
      }
      line=trim(line);
      lineStream = new IStringStream(line);
      if(currentElements_<0) { // Ordered type
        if (block_flag) {
          int max_points_ = currentI_*currentJ_*currentK_;
          // DATAPACKING = BLOCK:
          // Read data: x_0 ... x_max_points_ y_0 ... y_max_points ...
          scalar value_;
          while(*lineStream) {
            *lineStream >> value_;
            if (current_variable_>=data_.size()) {
              FatalErrorInFunction << "Problem reading the Tecplot file" << exit(FatalError);
            }
            data_[current_variable_].append(value_);
            current_points_number_++;
            if(current_points_number_>= max_points_) {
              current_variable_++;
              current_points_number_=0;
            }
          }
        }
        if (point_flag) {
          // DATAPACKING = POINT:
          // Read data: x_0 y0 ... x_max_points_ y_max_points ...
          scalar value_;
          int variables_size = variables_.size();
          while(*lineStream) {
            *lineStream >> value_;
            data_[current_variable_].append(value_);
            current_variable_++;
            if (current_variable_>=variables_size) {
              current_variable_=0;
            }
          }
        }
      } else { // FETriangle type
        // DATAPACKING = BLOCK:

        // Read nodes-elements: elem1_node1 elem1_node2 elem1_node3 elem2_node1 ...
        if (current_variable_>=data_.size()) {
          int value1_; // index first node for elem i
          int value2_; // index second node for elem i
          int value3_; // index third node for elem i
          while(*lineStream) {
            *lineStream >> value1_;
            *lineStream >> value2_;
            *lineStream >> value3_;
            vector v_(value1_-1,value2_-1, value3_-1); // elements indexes start at 1 in tecplot
            nodesElements_.append(v_);
          }
        } else {
          // Read data: x_0 ... x_max_points y_0 ... y_max_points ...
          scalar value_;
          int max_points_ = currentNodes_;  // READING nodes
          if (foundInList<scalar>(current_variable_,varlocation_)) {
            max_points_ = currentElements_; // READING elements
          }
          while(*lineStream) {
            *lineStream >> value_;
            if (current_variable_>=data_.size()) {
              break;
            }
            data_[current_variable_].append(value_);
            current_points_number_++;
            if(current_points_number_>= max_points_) {
              current_variable_++;
              current_points_number_=0;
            }
          }
        }
      }
    }
  }



  // Compute the centroid of the triangles using the nodes xyz
  if(nElements_>0) { // Ordered type
    List<vector> xyzNodes_(nNodes_);
    if (data_.size()<4) {
      FatalErrorInFunction << "Number of variables has to be more than 3. (x y z f1 f2 ...)" << exit(FatalError);
    }
    if ((data_[0].size()!= nNodes_)||(data_[1].size()!= nNodes_)||(data_[2].size()!= nNodes_)) {
      FatalErrorInFunction << "The 3 first variables must be X Y Z (nodes)." << exit(FatalError);
    }

    forAll(xyzNodes_, nodeI) {
      vector v_(data_[0][nodeI],data_[1][nodeI],data_[2][nodeI]);
      xyzNodes_[nodeI]=v_;
    }
    List<vector> centroids_ = centroidTriangle(nodesElements_,xyzNodes_);
    if(centroids_.size()!=nElements_) {
      FatalErrorInFunction << "Problem during the centroid calculation." << exit(FatalError);
    }
    data_[0].resize(nElements_);
    data_[1].resize(nElements_);
    data_[2].resize(nElements_);
    forAll(data_[0], elemI) {
      data_[0][elemI] = centroids_[elemI][0];
      data_[1][elemI] = centroids_[elemI][1];
      data_[2][elemI] = centroids_[elemI][2];
    }
  }

  if (data_.size()==0) {
    Info << filename_ << " is empty." << endl;
  } else {
    Info << "Read " << filename_ << " : data = (" << data_.size() << "x" << data_[0].size() << ")" << endl;
    Info << " variables = [ " ;
    forAll(variables_, varI) {
      if (varI==variables_.size()-1) {
        Info << varI << ":\"" << variables_[varI] << "\" ]";
      } else {
        Info << varI << ":\"" << variables_[varI] << "\" , ";
      }
    }
    Info << endl;
  }

  return data_;
}

////- Write Tecplot data file
//inline void writeTecplotFileData(const fileName filename_, const List<scalarList> data_, const word format_)
//{
//  if (format_ != "FETriangle" || format_ != "Ordered"){
//    FatalErrorInFunction << "format has to be \"FETriangle\" or \"Ordered\"" << exit(FatalError);
//  }
//  fileName fileName_ = changeEnviVar(filename_);
//  OFstream dataFile(fileName_);
//  if (!dataFile.good()) {
//    FatalErrorInFunction
//        << "Cannot write file " << filename_
//        << exit(FatalError);
//  }
//  // Writing data file
//  dataFile << "TITLE = \"title\"" << endl;
//  dataFile << "VARIABLES = " << endl;
//  forAll(data, i){
//    dataFile << "\"var-" << i << "\"" << endl;
//  }
//  dataFile << "ZONE T=\"ZONE 001\"" << endl;
//
//
//}


//- Gives the list of sub folders path with a given keyword
inline void searchSubFoldersKeyword(List<fileName>& listFileName, fileName sourceFolder, word keyword)
{
  sourceFolder=changeEnviVar(sourceFolder);
  fileNameList dirs=readDir(sourceFolder, fileType::directory);
  forAll(dirs, i) {
    fileName path = search(keyword, sourceFolder/dirs[i]);
    if (path!="") {
      bool addPath=true;
      forAll(listFileName, pathI) {
        if(path.path()==listFileName[pathI]) {
          addPath=false;
        }
      }
      if (addPath) {
        listFileName.append(path.path());
      }
    }
    searchSubFoldersKeyword(listFileName, sourceFolder/dirs[i],keyword);
  }
}

//- Gives the list of folders path with a given keyword
inline List<fileName> searchFoldersKeyword(fileName sourceFolder, word keyword)
{
  sourceFolder=changeEnviVar(sourceFolder);
  List<fileName> foldersKeyword_;
  if (!isDir(sourceFolder)) {
    FatalErrorInFunction << sourceFolder << " not found." << exit(FatalError);
  }
  fileName path = search(keyword, sourceFolder);
  if (path!="") {
    foldersKeyword_.append(path.path());
  }
  searchSubFoldersKeyword(foldersKeyword_,sourceFolder,keyword);
  return foldersKeyword_;
}


inline void filesInSubFolder(List<fileName>& listFileName, fileName sourceFolder)
{
  sourceFolder=changeEnviVar(sourceFolder);
  fileNameList dirs=readDir(sourceFolder, fileType::directory);
  fileNameList files=readDir(sourceFolder, fileType::file);
  forAll(files, i) {
    bool addPath=true;
    forAll(listFileName, pathI) {
      if(files[i]==listFileName[pathI]) {
        addPath=false;
      }
    }
    if (addPath) {
      listFileName.append(sourceFolder/files[i]);
    }
  }

  forAll(dirs, i) {
    filesInSubFolder(listFileName, sourceFolder/dirs[i]);
  }

}


//- Gives the list of all the files in a folder
inline List<fileName> filesInFolder(fileName sourceFolder)
{
  sourceFolder=changeEnviVar(sourceFolder);
  if (!isDir(sourceFolder)) {
    FatalErrorInFunction << sourceFolder << " not found." << exit(FatalError);
  }
  List<fileName> listFileName;
  filesInSubFolder(listFileName, sourceFolder);
  return listFileName;
}

inline void printDescription( const fvMesh& mesh, const word type)
{
  IOdictionary dict
  (    IOobject
       (
           "fieldsDescription",
           (word) changeEnviVar("$PATO_DIR")+"/documentation/Fields",
           mesh,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       )
  );
  const wordList list(mesh.objectRegistry::sortedNames(type));
  if (dict.isDict(type)) {
    dictionary subdict(dict.subDict(type));
    FatalErrorInFunction << list.size() << nl << "(" << endl;
    forAll(list, i) {
      string print(subdict.lookupOrDefault<string>(list[i], "Description not found in $PATO_DIR/documentation/Fields/fieldsDescription"));
      FatalErrorInFunction << list[i] << " = " << print << endl;
    }
    FatalErrorInFunction << ")" << endl;
  } else {
    FatalErrorInFunction << list << endl;
  }
}

//- Returns true if the fields are found in mesh
inline void foundFieldsInMesh
(
    const fvMesh& mesh_,
    const wordList nameFieldsList_
)
{
//  wordList newList_;
  int indexNotFound=-1;
  forAll(nameFieldsList_, wordI) {
    if (!mesh_.objectRegistry::foundObject<volScalarField>(nameFieldsList_[wordI])) {
      if (!mesh_.objectRegistry::foundObject<volVectorField>(nameFieldsList_[wordI])) {
        if (!mesh_.objectRegistry::foundObject<volTensorField>(nameFieldsList_[wordI])) {
          if (!mesh_.objectRegistry::foundObject<volSymmTensorField>(nameFieldsList_[wordI])) {
            if(!mesh_.objectRegistry::foundObject<surfaceScalarField>(nameFieldsList_[wordI])) {
              if(!mesh_.objectRegistry::foundObject<surfaceVectorField>(nameFieldsList_[wordI])) {
                if(!mesh_.objectRegistry::foundObject<surfaceTensorField>(nameFieldsList_[wordI])) {
                  indexNotFound = wordI;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  if (foundInList<word>("all",nameFieldsList_)) {
    FatalErrorInFunction << "\"all\" fields can not be used." << exit(FatalError);
  }
  if(indexNotFound>=0) {
    FatalErrorInFunction
        << nameFieldsList_[indexNotFound] << " not found in the mesh fields." << nl
        << "\e[1mvolScalarField available: \e[0m ";
    printDescription(mesh_, "volScalarField");
    FatalErrorInFunction << nl <<  "\e[1msurfaceScalarField available: \e[0m";
    printDescription(mesh_, "surfaceScalarField");
    FatalErrorInFunction << nl <<  "\e[1mvolVectorField available: \e[0m";
    printDescription(mesh_, "volVectorField");
    FatalErrorInFunction << nl <<  "\e[1msurfaceVectorField available: \e[0m";
    printDescription(mesh_, "surfaceVectorField");
    FatalErrorInFunction << nl <<  "\e[1mvolTensorField available: \e[0m";
    printDescription(mesh_, "volTensorField");
    FatalErrorInFunction << nl <<  "\e[1msurfaceTensorField available: \e[0m";
    printDescription(mesh_, "surfaceTensorField");
    FatalErrorInFunction << nl <<  "\e[1mvolSymmTensorField available: \e[0m";
    printDescription(mesh_, "volSymmTensorField");
    FatalErrorInFunction << exit(FatalError);
  }
}

//- Returns true if the fields are found in mesh
inline void foundFieldsInMeshAll
(
    const fvMesh& mesh_,
    wordList& nameFieldsList_
)
{
//  wordList newList_;
  int indexNotFound=-1;
  forAll(nameFieldsList_, wordI) {
    if (!mesh_.objectRegistry::foundObject<volScalarField>(nameFieldsList_[wordI])) {
      if (!mesh_.objectRegistry::foundObject<volVectorField>(nameFieldsList_[wordI])) {
        if (!mesh_.objectRegistry::foundObject<volTensorField>(nameFieldsList_[wordI])) {
          if (!mesh_.objectRegistry::foundObject<volSymmTensorField>(nameFieldsList_[wordI])) {
            if(!mesh_.objectRegistry::foundObject<surfaceScalarField>(nameFieldsList_[wordI])) {
              if(!mesh_.objectRegistry::foundObject<surfaceVectorField>(nameFieldsList_[wordI])) {
                if(!mesh_.objectRegistry::foundObject<surfaceTensorField>(nameFieldsList_[wordI])) {
                  indexNotFound = wordI;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  if (foundInList<word>("all",nameFieldsList_)) {
    indexNotFound=-1;
    nameFieldsList_.clear();
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("volScalarField"));
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("volVectorField"));
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("volTensorField"));
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("volSymmTensorField"));
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("surfaceScalarField"));
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("surfaceVectorField"));
    nameFieldsList_.append(mesh_.objectRegistry::sortedNames("surfaceTensorField"));
    wordList cleanedAllFields;
    forAll(nameFieldsList_, fieldI) {
      if(nameFieldsList_[fieldI]!="") {
        cleanedAllFields.append(nameFieldsList_[fieldI]);
      }
    }
    nameFieldsList_=cleanedAllFields;
  }
  if(indexNotFound>=0) {
    FatalErrorInFunction
        << nameFieldsList_[indexNotFound] << " not found in the mesh fields." << nl
        << "\e[1mvolScalarField available: \e[0m ";
    printDescription(mesh_, "volScalarField");
    FatalErrorInFunction << nl <<  "\e[1msurfaceScalarField available: \e[0m";
    printDescription(mesh_, "surfaceScalarField");
    FatalErrorInFunction << nl <<  "\e[1mvolVectorField available: \e[0m";
    printDescription(mesh_, "volVectorField");
    FatalErrorInFunction << nl <<  "\e[1msurfaceVectorField available: \e[0m";
    printDescription(mesh_, "surfaceVectorField");
    FatalErrorInFunction << nl <<  "\e[1mvolTensorField available: \e[0m";
    printDescription(mesh_, "volTensorField");
    FatalErrorInFunction << nl <<  "\e[1msurfaceTensorField available: \e[0m";
    printDescription(mesh_, "surfaceTensorField");
    FatalErrorInFunction << nl <<  "\e[1mvolSymmTensorField available: \e[0m";
    printDescription(mesh_, "volSymmTensorField");
    FatalErrorInFunction << exit(FatalError);
  }
}

// number of columns
inline int numberColumns(const string& line)
{

  string::size_type pos = 0;
  int nNumbers = 0;
  string s = line;
  string numbers;
  string delimiter=" ";
  while ((pos = s.find(delimiter)) != std::string::npos) {
    numbers = s.substr(0, pos);

    if (numbers != "") {
      nNumbers++;
    }
    s.erase(0, pos + delimiter.length());
  }
  if (s != "" ) {
    nNumbers++;
  }
  return nNumbers;
}

// number of scalar columns
inline int numberColumnsScalar(string& line)
{
  string::size_type pos = 0;
  int nNumbers = 0;
  string s = line;
  string numbers;
  string delimiter=" ";
  while ((pos = s.find(delimiter)) != std::string::npos) {
    numbers = s.substr(0, pos);

    if (numbers != "") {
      if (isNumber(numbers)) {
        nNumbers++;
      }
    }
    s.erase(0, pos + delimiter.length());
  }
  if (s != "" ) {
    if (isNumber(s)) {
      nNumbers++;
    }
  }
  return nNumbers;
}

template <typename Type> Type& meshLookupOrConstructModel
(
    const fvMesh& mesh,
    const word& regionName,
    const word& modelName
)
{
  const IOdictionary materialDict(
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
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  word modelModelName= modelName+"Model";
  word modelTypeName= "no";
  if (materialDict.isDict(modelName)) {
    dictionary modelDict = materialDict.subDict(modelName );
    modelTypeName = word(modelDict.lookupOrDefault<word>(modelName+"Type","no"));
  }
  if (!mesh.objectRegistry::foundObject<Type>(regionName+modelName+"Model")) {
    // Print the message
    Info << Type::getTabLevel("| ") << bold_on << "Create "<< (string) modelModelName
         << " type " << (string) modelTypeName << bold_off << endl;
    // Create the model
    autoPtr<Type> fPtr(Type::New(mesh, regionName));
    fPtr->store(fPtr);
  } else {
    if (Type::tabLevel_<=0) {
      Type::tabLevel_=0;
      Info << Type::getTabLevel("| ");
    } else {
      Info << Type::getTabLevel();
    }
    Info << "Reference to " << (string) modelModelName
         << " type " << (string) modelTypeName << endl;
  }

  return const_cast<Type&>
         (
             mesh.objectRegistry::lookupObject<Type>(regionName+modelName+"Model")
         );
}

template<typename Type>
inline UniformDimensionedField<Type>& meshUniqueConstruct
(
    const fvMesh& mesh,
    const word& name,
    const word modelName = ""
)
{
  if (mesh.objectRegistry::foundObject<UniformDimensionedField<Type>>(name)) {
    if (modelName!="") {
      FatalError << modelName << "Model: ";
    }
    FatalError << (fileName) name << " exists already in mesh " << (fileName) mesh.name() << exit(FatalError);
  } else {
    UniformDimensionedField<Type>* fPtr
    (
        new UniformDimensionedField<Type>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<UniformDimensionedField<Type>&>
         (
             mesh.objectRegistry::lookupObject<UniformDimensionedField<Type>>(name)
         );
}

template<typename Type1, template<typename> class Type2, class Type3>
inline GeometricField<Type1, Type2, Type3>& meshUniqueConstruct
(
    const fvMesh& mesh,
    const word& name,
    const word modelName = ""
)
{
  if (mesh.objectRegistry::foundObject<GeometricField<Type1, Type2, Type3>>(name)) {
    if (modelName!="") {
      FatalError << modelName << "Model: ";
    }
    FatalError << (fileName) name << " exists already in mesh " << (fileName) mesh.name() << exit(FatalError);
  } else {
    GeometricField<Type1, Type2, Type3>* fPtr
    (
        new GeometricField<Type1, Type2, Type3>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<GeometricField<Type1, Type2, Type3>&>
         (
             mesh.objectRegistry::lookupObject<GeometricField<Type1, Type2, Type3>>(name)
         );
}

template<typename Type>
inline UniformDimensionedField<Type>& meshUniqueConstruct
(
    const fvMesh& mesh,
    const word& name,
    const dimensioned<Type>& dim_,
    const word modelName = "",
    const IOobject::readOption& readOpt_ = IOobject::READ_IF_PRESENT
)
{
  if (mesh.objectRegistry::foundObject<UniformDimensionedField<Type>>(name)) {
    if (modelName!="") {
      FatalError << modelName << "Model: ";
    }
    FatalError << (fileName) name << " exists already in mesh " << (fileName) mesh.name() << exit(FatalError);
  } else {
    UniformDimensionedField<Type>* fPtr
    (
        new UniformDimensionedField<Type>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                readOpt_,
                IOobject::NO_WRITE
            ),
            dim_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<UniformDimensionedField<Type>&>
         (
             mesh.objectRegistry::lookupObject<UniformDimensionedField<Type>>(name)
         );
}

template<typename Type1, template<typename> class Type2, typename Type3>
inline GeometricField<Type1, Type2, Type3>& meshUniqueConstruct
(
    const fvMesh& mesh,
    const word& name,
    const dimensioned<Type1>& dim_,
    const word modelName = "",
    const IOobject::readOption& readOpt_ = IOobject::READ_IF_PRESENT
)
{
  if (mesh.objectRegistry::foundObject<GeometricField<Type1, Type2, Type3>>(name)) {
    if (modelName!="") {
      FatalError << modelName << "Model: ";
    }
    FatalError << (fileName) name << " exists already in mesh " << (fileName) mesh.name() << exit(FatalError);
  } else {
    GeometricField<Type1, Type2, Type3>* fPtr
    (
        new GeometricField<Type1, Type2, Type3>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                readOpt_,
                IOobject::NO_WRITE
            ),
            mesh,
            dim_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<GeometricField<Type1, Type2, Type3>&>
         (
             mesh.objectRegistry::lookupObject<GeometricField<Type1, Type2, Type3>>(name)
         );
}

template<typename Type1, template<typename> class Type2, typename Type3>
inline GeometricField<Type1, Type2, Type3>& meshUniqueConstruct
(
    const fvMesh& mesh,
    const word& name,
    const dimensioned<Type1>& dim_,
    const wordList& BC_,
    const word modelName = "",
    const IOobject::readOption& readOpt_ = IOobject::READ_IF_PRESENT
)
{
  if (mesh.objectRegistry::foundObject<GeometricField<Type1, Type2, Type3>>(name)) {
    if (modelName!="") {
      FatalError << modelName << "Model: ";
    }
    FatalError << (fileName) name << " exists already in mesh " << (fileName) mesh.name() << exit(FatalError);
  } else {
    GeometricField<Type1, Type2, Type3>* fPtr
    (
        new GeometricField<Type1, Type2, Type3>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                readOpt_,
                IOobject::NO_WRITE
            ),
            mesh,
            dim_,
            BC_
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<GeometricField<Type1, Type2, Type3>&>
         (
             mesh.objectRegistry::lookupObject<GeometricField<Type1, Type2, Type3>>(name)
         );
}


template<typename Type1, template<typename> class Type2, typename Type3>
inline GeometricField<Type1, Type2, Type3>& meshUniqueConstruct
(
    const fvMesh& mesh,
    const word& name,
    const GeometricField<Type1, Type2, Type3>& field,
    const word modelName = "",
    const IOobject::readOption& readOpt_ = IOobject::READ_IF_PRESENT
)
{
  if (mesh.objectRegistry::foundObject<GeometricField<Type1, Type2, Type3>>(name)) {
    if (modelName!="") {
      FatalError << modelName << "Model: ";
    }
    FatalError << (fileName) name << " exists already in mesh " << (fileName) mesh.name() << exit(FatalError);
  } else {
    GeometricField<Type1, Type2, Type3>* fPtr
    (
        new GeometricField<Type1, Type2, Type3>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                readOpt_,
                IOobject::NO_WRITE
            ),
            field
        )
    );

    // Transfer ownership of this object to the objectRegistry
    fPtr->store(fPtr);
  }

  return const_cast<GeometricField<Type1, Type2, Type3>&>
         (
             mesh.objectRegistry::lookupObject<GeometricField<Type1, Type2, Type3>>(name)
         );
}

template<typename Type1, template<typename> class Type2, typename Type3>
inline GeometricField<Type1, Type2, Type3>& meshLookup
(
    const fvMesh& mesh,
    const word name,
    const word modelName=""
)
{
  if (!mesh.objectRegistry::foundObject<GeometricField<Type1, Type2, Type3>>(name)) {
    FatalError << name << " not found in mesh " << mesh.name() << ".";
    if (modelName!="") {
      FatalError << "Please choose another " << modelName << "Type.";
    }
    FatalError << exit(FatalError);
  }

  return const_cast<GeometricField<Type1, Type2, Type3>&>
         (
             mesh.objectRegistry::lookupObject<GeometricField<Type1, Type2, Type3>>(name)
         );
}

template<typename Type>
inline Type& meshLookup
(
    const fvMesh& mesh,
    const word name,
    const word modelName
)
{
  if (!mesh.objectRegistry::foundObject<Type>(name)) {
    FatalError << name << " not found in mesh " << mesh.name() << ".";
    if (modelName!="") {
      FatalError << "Please choose another " << modelName << "Type.";
    }
    FatalError << exit(FatalError);
  }

  return const_cast<Type&>
         (
             mesh.objectRegistry::lookupObject<Type>(name)
         );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
