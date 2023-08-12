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
  patoCodingStyle

  Description
  Format the code using google-astyle.

  \*---------------------------------------------------------------------------*/

// include headers
#include "fvCFD.H"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>


// namespace
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return the executation of a command
std::string exec(const char* cmd)
{
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

// Return all the files in the directory and sub-directories
fileNameList files_in_dir(fileName dirPath, fileNameList& all_files, Switch verbose=false)
{
  if (isDir(dirPath)) {
    if (dirPath.name()=="Make" || dirPath.name()=="lnInclude") {
      return fileNameList();
    }
  }
  const fileName pato_dir=getEnv("PATO_DIR");
  fileNameList hidden_dir_pato= {
    pato_dir/"documentation",
    pato_dir/"src/thirdParty/mutation++",
    pato_dir/"install"
  };
  forAll(hidden_dir_pato, i) {
    if (dirPath==hidden_dir_pato[i]) {
      return fileNameList();
    }
  }
  if (verbose)
    Info << "readDir(" << dirPath << ", fileType::file);" << endl;
  fileNameList files=readDir(dirPath, fileType::file, false);
  forAll(files, i) {
    if (files[i]!="README.md")
      all_files.append(dirPath/files[i]);
  }
  if (verbose)
    Info << "readDir(" << dirPath << ", fileType::directory);" << endl;
  fileNameList dirs=readDir(dirPath, fileType::directory);
  forAll(dirs, i) {
    files_in_dir(dirPath/dirs[i], all_files, verbose);
  }
  return files;
}

// Check if a file is valid for the formatting (C++ files only)
bool valid_file(fileName file)
{
  // Check extensions
  if (file.ext()=="h" || file.ext()=="H" || file.ext()=="C")
    return true;
  if (file.ext() == "gz" || file.ext() == "PNG" || file.ext() == "xls")
    return false;

  // Check if the first not empty line of the file starts with the "#" character
  word cmd = "echo \"`awk '/./{print;exit}' \"" + file + "\"`\" | head -c 1";
  word s = exec(cmd.c_str());
  if (s == "#") {
    return false;
  }

  // Check if there is an OpenFOAM header
  Foam::Time runTime(file.path(),file.name());
  IOobject io
  (
      "",
      "",
      runTime.db(),
      IOobject::MUST_READ,
      IOobject::NO_WRITE
  );
  bool ok = true;
  int level=Warning.level;
  Warning.level=0;
  ok = Foam::fileHandler().readHeader(io, file, "");
  Warning.level=level;
  if (ok)
    return true;

  return false;
}

int main(int argc, char *argv[])
{
  // Options
  argList::noParallel();
  argList::removeOption("case");
  argList::removeOption("noFunctionObjects");
  argList::addOption("verbose","","Print all the files.");
  argList::addOption("format","Switch","Format the code. Default='true'.");
  argList::addOption("dir","path","Path of the directory where the files will be formatted. Default='cwd()'");

  argList args(argc, argv);

  if (!args.check()) {
    FatalError.exit();
  }

  Switch verbose = false;
  if(args.optionReadIfPresent("verbose", verbose)) {
    verbose = true;
  }

  Switch format = true;
  args.optionReadIfPresent("format", format);

  fileName dir_path = getEnv("PWD");
  args.optionReadIfPresent("dir", dir_path);

  Info << nl << "Reading " << dir_path << nl << endl;
  int n_total=0;
  int n_formatted=0;
  fileNameList files;
  if (verbose)
    Info << "files_in_dir(dir_path, files, verbose):" << endl;
  files_in_dir(dir_path, files, verbose);
  if (verbose)
    Info << nl << "files = " << files << endl;
  fileNameList valid_files;

  forAll(files, i) {
    if (valid_file(files[i])) {
      valid_files.append(files[i]);
      if (verbose && !format) {
        word start = "  $ ";
        if (!format)
          start += "#";
        Info << start << "style-google " << files[i] << " --style=linux --indent=spaces=2" << endl;
      }
      n_total++;
      if (format) {
        word cmd = ". $PATO_DIR/src/applications/utilities/runFunctions/RunFunctions; pato_init; style-google \""
                   +files[i]+"\" --style=linux --indent=spaces=2";
        word s = exec(cmd.c_str());
        if (s.substr(0,9) == "Formatted") {
          n_formatted++;
          if (n_formatted==1) {
            Info << nl << "Formatted files:" << endl;
          }
          if (isFile(files[i]+".orig")) {
            word cmd = "rm \""+files[i]+"\".orig";
            system(cmd);
          }
          Info << "  - " << files[i] << endl;
        }
      }
    }
  }

  if (verbose && n_formatted!=0)
    Info << endl;
  if (verbose)
    Info << "valid files = " << valid_files << endl;

  std::string bold = "\e[1m";
  std::string regular = "\e[0m";
  if (n_formatted>0) {
    Info << endl;
    std::cout << "Total number of valid files = " << n_total << std::endl;
    std::cout << bold << "Total number of formatted files = "
              << n_formatted << regular << std::endl;
    Info << endl;
  } else {

    std::cout << "Total number of valid files = " << n_total << std::endl;
    std::cout << bold << "No formatted files" << regular << std::endl;
    Info << endl;
  }

  return 0;
}
