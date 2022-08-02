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
  along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  Utility
  postProcessTime

  Description
  Create a file (field values in function of time) from the postProcessing files.  

  \*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PATOx.H" // PATOx headers
#include "IOFunctions.H" // IO functions headers
#include "IOmanip.H" // IO manip headers
#include "fileOperation.H" // file operation headers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
  argList::noParallel();
  argList::validArgs.append("dictName");
  argList::validArgs.append("regionName");
  argList::validArgs.append("inputFilename");
  argList::validArgs.append("outputFilename");
  argList args(argc, argv);

  if (!args.check()) {
    FatalError.exit();
  }

  word dictName_ = args[1]; // dict name
  word regionName_ = args[2]; // region name 
  fileName file_ = args[3]; // input file name    
  fileName outputName_= args[4]; // output file name
  OFstream os_(outputName_); // output stream
  fileName dirName_ = "postProcessing/"+dictName_+"/"+regionName_;
  
  // Verify dirName_ exists
  if (!isDir(dirName_)){
    FatalError << dirName_ << " not found." << exit(FatalError);
  }

  // Name of the time directories
  fileNameList dirTimes=readDir(dirName_, fileType::directory);
  scalarList times;

  forAll(dirTimes, dirI){
    scalar time = std::stof(dirTimes[dirI]);
    times.append(time);
  }

  // Put the times in sorted order
  labelList visitOrder;
  sortedOrder(times, visitOrder);
  times = scalarList(times, visitOrder);
  dirTimes = fileNameList(dirTimes, visitOrder);

  // Write the time data in output file
  forAll(dirTimes, dirI) {
    fileName fileI_ = dirName_ +"/"+dirTimes[dirI]+ "/"+ file_.name();
    List<scalarList> dataI_(readFileData(fileI_));
    if (dirI==0){
      const int n_fields = dataI_.size()-3;
      if (n_fields <= 0){
	FatalError << "No fields found in " << fileI_ << exit(FatalError);
      }
      os_<< "//t(s) ";
      for (int i = 0; i < n_fields; i++){
	os_ << "field_"<<i+1<< " ";
      }
      os_ << endl;
    }
    os_<< (word) dirTimes[dirI] << " ";
    forAll(dataI_, fieldI) {
      if (fieldI >2) { // X Y Z
        forAll(dataI_[fieldI], probI) {
          os_<< dataI_[fieldI][probI] << " ";
        }
      }
    }
    os_ << endl;
  }
  Info << "Writing " << outputName_ << endl;


  return 0;
}


// ************************************************************************* //
