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

  Application
  replaceDakotaResults

  Description
  Replace the PATO template file using the Dakota results

  \*---------------------------------------------------------------------------*/

// include headers
#include "fvCFD.H"
#include "unitConversionBritish.H"
#include "mathFunctions.H"

// namespace
using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
  argList::noParallel();
  argList::validArgs.append("results");
  argList::validArgs.append("template");
  argList::validArgs.append("fileToChange");

  argList args(argc, argv);


  if (!args.check()) {
    FatalError.exit();
  }

  word resultsFileName = args.args()[1];

  Info << "Reading " << resultsFileName << endl;

  IFstream dataFile(resultsFileName);
  PtrList<OFstream> outputFiles;

  if (!dataFile.good()) {
    FatalErrorInFunction
        << "Cannot read file " << resultsFileName
        << exit(FatalError);
  }


  int lineI = 0;

  wordList titles_;
  scalarList results_;
  IStringStream * lineStream;

  while(dataFile) {
    word line;
    dataFile.getLine(line);
    if (line.empty()) { // empty line or only space/tab
      break;
    }

    lineStream = new IStringStream(line);
    if (lineI == 2) {
      word a = (word) *lineStream; // " "
      a = (word) *lineStream; // "Row"
      a = (word) *lineStream; // "Labels:"
      while(!a.empty()) {
        titles_.append(a);
        a = (word) *lineStream;
      }
    }
    if(lineI>4 && lineI <= titles_.size()+4) {
      results_.append(readScalar(*lineStream));
    }

    lineI++;
  }

  // sed command
  word sed_cmd = "sed";
#if defined(darwin64)
  sed_cmd = "gsed";
#endif
  fileName templateFile_ = args.args()[2];
  fileName fileToChange_ = args.args()[3];
  if(isFile(templateFile_)) {
    word command = "cp " + templateFile_ + " " + fileToChange_;
    Info << "system: " << command << endl;
    system(command);
    forAll(titles_, titleI) {
      command = sed_cmd + " -i \"s/{" + titles_[titleI] + "}/" + name(results_[titleI]) + "/g\"" + " " + fileToChange_ ;
      Info << "system: " << command << endl;
      system(command);
    }
  } else {
    FatalErrorInFunction << templateFile_ << " not found." << exit(FatalError);
  }

  return 0;
}
