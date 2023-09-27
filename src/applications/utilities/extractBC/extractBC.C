/*---------------------------------------------------------------------------* \
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
  extractBC

  Description
  Extract the BC to a txt file (X Y Z field)

  \*---------------------------------------------------------------------------*/

// Macro for Info
#define Info_name(name) Info << #name << "=\"" << name << "\"" << endl

#include "fvCFD.H"
#include "IOmanip.H"
#include <filesystem>

// Return all the files in the directory and sub-directories
fileNameList files_in_dir(fileName dirPath, fileNameList& all_files)
{
  fileNameList files=readDir(dirPath, fileType::file);
  forAll(files, i) {
    all_files.append(dirPath/files[i]);
  }
  fileNameList dirs=readDir(dirPath, fileType::directory);
  forAll(dirs, i) {
    files_in_dir(dirPath/dirs[i], all_files);
  }
  return files;
}

// Return the files with an extension
fileNameList files_ext(const fileNameList& all_files, word ext)
{
  fileNameList new_files;
  forAll(all_files, fileI) {
    if(all_files[fileI].name().substr(all_files[fileI].name().find_last_of(".")) == ext) {
      new_files.append(all_files[fileI]);
    }
  }
  return new_files;
}

// Return the files with a substring
fileNameList files_substr(const fileNameList& all_files, word substring)
{
  fileNameList new_files;
  forAll(all_files, fileI) {
    int found = all_files[fileI].name().find(substring);
    if(found>=0) {
      new_files.append(all_files[fileI]);
    }
  }
  return new_files;
}

// Return all the lines of a file
wordList read_lines(fileName path)
{
  // Reading data file
  wordList lines;
  IFstream dataFile(path);
  if (!dataFile.good()) {
    FatalErrorInFunction
        << "Cannot read file " << path
        << exit(FatalError);
  }
  while(dataFile) {
    string line;
    dataFile.getLine(line);
    lines.append(line);
  }
  return lines;
}

// Extract a string between delimiters in a file
word read_str_delim(fileName path, word delim_left, word delim_right)
{
  word str="";
  wordList lines = read_lines(path);
  forAll(lines, lineI) {
    string line = lines[lineI];
    // remove space
    int len = line.length(); // storing the length of the string
    int count = std::count(line.begin(), line.end(), ' '); // counting the number of whitespaces
    std::remove(line.begin(), line.end(),' '); // removing all the whitespaces
    line.resize(len - count); // resizing the string to len-count
    int found_left=line.find(delim_left);
    int found_right=line.find(delim_right);
    if (found_left >=0 && found_right >=0) {
      str=line.substr(found_left+delim_left.size(),found_right-(found_left+delim_left.size()));
      break;
    }
  }
  return str;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
  // Executable arguments
  argList::noParallel();
  argList::validArgs.append("inputFile");
  argList args(argc, argv);
  if (!args.check()) {
    FatalError.exit();
  }

  // Input file dictionary
  fileName file_ = args[1];
  Foam::Time runTime2(file_.path(),file_.name());
  IOdictionary dict(
      IOobject
      (
          "",
          "",
          runTime2.db(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  // Inputs
  const word regionName=dict.lookup("regionName"); // Region name
  const word patchName=dict.lookup("patchName"); // Patch name
  const word fieldName=dict.lookup("fieldName"); // Field name
  const word units=dict.lookup("units"); // Field units [...]
  const scalarList times=dict.lookup("times"); // Time [s]
  const word outputDir=dict.lookup("outputDir"); // Output name

  // Check output dir
  if (isDir(outputDir)) {
    FatalError << outputDir << " exists already." << exit(FatalError);
  }
  Info << "Create " << outputDir << endl;
  system("mkdir -p "+outputDir);

  // Create time/field and write output dir
#include "createTime.H"
  instantList timeDirs = timeSelector::select0(runTime, args);
  forAll(times, timeI) {
    // Check if time exists
    scalar time = times[timeI];
    int index_t=-1;
    forAll(timeDirs, timeJ) {
      if (timeDirs[timeJ].value()==time) {
        index_t=timeJ;
        break;
      }
    }
    if (index_t<0) {
      FatalError << "time (" << time << ") not found." << exit(FatalError);
    }

    // Get all the boundary names in PATO
    const fileName PATO_DIR = getEnv("PATO_DIR");
    const fileName BC_path = PATO_DIR+"/src/applications/libraries/libPATOx/MaterialModel/BoundaryConditions";
    fileNameList files;
    files_in_dir(BC_path, files);
    files = files_ext(files,".H");
    files = files_substr(files,"FvPatch");
    wordList BC_names;
    forAll(files, fileI) {
      BC_names.append(read_str_delim(files[fileI],"TypeName(\"","\");"));
    }

    // Replace the boundary names from PATO to OpenFOAM
    const fileName field_path = getEnv("FOAM_CASE")+"/"+name(times[timeI])+"/"+regionName+"/"+fieldName;
    Info << "Replace the PATO BC to fixedValue in " << field_path << endl;
    word sed_cmd = "sed"; // sed command
#if defined(darwin64)
    sed_cmd = "gsed";
#endif
    forAll(BC_names, nameI) {
      word cmd = sed_cmd + " -i \"s/type \\+"+ BC_names[nameI]+" \\+;/type fixedValue;/g\" "+field_path;
      system(cmd);
      cmd = sed_cmd + " -i \"s/type \\+"+ BC_names[nameI]+";/type fixedValue;/g\" "+field_path;
      system(cmd);
    }

    // Change time
    runTime.setTime(timeDirs[index_t], index_t);

    // Create mesh
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    // Get patchID
    const label patchID=mesh.boundaryMesh().findPatchID(patchName);

    // Create field
    volScalarField field(
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Write output
    fileName outputPath = outputDir+"/output_"+name(time);
    Info << "Writing " << outputPath << nl << endl;
    OFstream os_out(outputPath);
    os_out.setf(ios_base::scientific, ios_base::floatfield);
    os_out.setf(ios_base::left);
    os_out.precision(5);
    wordList header= {"// x[m]","y[m]","z[m]",fieldName+units};
    forAll(header,i) {
      os_out << setw(16) << header[i];
    }
    os_out << endl;
    forAll(mesh.Cf().boundaryField()[patchID], faceI) {
      os_out << setprecision(6) << setw(16) << mesh.Cf().boundaryField()[patchID][faceI][0];
      os_out << setprecision(6) << setw(16) << mesh.Cf().boundaryField()[patchID][faceI][1];
      os_out << setprecision(6) << setw(16) << mesh.Cf().boundaryField()[patchID][faceI][2];
      os_out << setprecision(6) << setw(16) << field.boundaryField()[patchID][faceI];
      os_out << endl;
    }
  }
  return 0;
}

// ************************************************************************* //
