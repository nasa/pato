#include "subtest.h"
#include <vector>
#include "PATOx.H"

class TestRestart : public SubTest
{
 public:

  TestRestart() {
    testSuiteName = "TestRestart";
    tests.push_back(test1);
    tests.push_back(test2);
  }

  // Run test1: comparison between AblationTestCase_1.x origin and restart cases
  static TestResult test1() {
    return run_test("test1");
  }

  // Run test1: comparison between AblationTestCase_2.x origin and restart cases
  static TestResult test2() {
    return run_test("test2");
  }

  // Run test: comparison between origin and restart cases
  static TestResult run_test(string test_name) {
    // Mute Info
    int level_Info_ = Info.level;
    Info.level = 0;
    // Change the tolerance for Linux
#ifndef __APPLE__
    scalarList old_tol=SubTest::tol;
    scalar new_tol=1e-2;
    forAll(SubTest::tol, tolI) {
      SubTest::tol[tolI]=new_tol;
    }
#endif
    // Create TestResult
    std::string testName = "Restart: regression testing";
    std::string testDescription = "Run test1 and compare the restart results to the origin results";
    fileName pato_unit_testing = getEnv("PATO_UNIT_TESTING");
    fileName testFolder = pato_unit_testing+"/testsuites/global/restart/"+test_name;
    fileName originTestFolder = testFolder+"/origin";
    fileName restartTestFolder = testFolder+"/restart";
    std::string suiteName = "Restart";
    TestResult result(suiteName, testName, 1, testDescription);
    // Clean and run test folder
    fileName log_file = testFolder+"/log.run";
    system(testFolder+"/Allclean > "+log_file+" 2>&1");
    system(testFolder+"/Allrun > "+log_file+" 2>&1");
    // Get time folders
    fileNameList originDirs=readDir(testFolder+"/origin", fileType::directory);
    fileNameList restartDirs=readDir(testFolder+"/restart", fileType::directory);
    fileNameList originTimeDirs;
    fileNameList restartTimeDirs;
    forAll(originDirs, dirI) {
      if(isNumber(originDirs[dirI])) {
        originTimeDirs.append(originDirs[dirI]);
      }
    }
    forAll(restartDirs, dirI) {
      if(isNumber(restartDirs[dirI])) {
        restartTimeDirs.append(restartDirs[dirI]);
      }
    }
    // Sort the list of folders
    labelList visitOrder;
    sortedOrder(originTimeDirs, visitOrder);
    originTimeDirs = fileNameList(originTimeDirs, visitOrder);
    sortedOrder(restartTimeDirs, visitOrder);
    restartTimeDirs = fileNameList(restartTimeDirs, visitOrder);
    // Compare all the files in the common time folders
    word regionName="porousMat";
    int index_restart=0;
    forAll(originTimeDirs, dirI) {
      Info << "t=" << originTimeDirs[dirI] << " [s]" << endl;
      scalar originTime = stof(originTimeDirs[dirI]);
      if (index_restart>restartTimeDirs.size()-1) {
        FatalError << "index_restart>restartTimeDirs.size()-1" << nl
                   << "restartTimeDirs = " << restartTimeDirs << nl
                   << "index_restart = " << index_restart << exit(FatalError);
      }
      scalar restartTime = stof(restartTimeDirs[index_restart]);
      if (originTime < restartTime) {
        continue;
      }
      // Check time folders
      if(!assertEquals(originTime, restartTime, &result)) {
        result.expected="originTime = " + name(originTime);
        result.expected+="\n";
        result.actual="restartTime = " + name(restartTime);
        result.actual+="\n";
        // Clean folder
        system(testFolder+"/Allclean > "+log_file+" 2>&1");
        system("rm -f "+testFolder+"/log.run");
        return result;
      }
      // Get all the sorted files in time folder
      fileNameList originFiles=filesInFolder(testFolder+"/origin/"+originTimeDirs[dirI]+"/"+regionName);
      fileNameList restartFiles=filesInFolder(testFolder+"/restart/"+restartTimeDirs[index_restart]+"/"+regionName);
      // Sort the list of folders
      labelList filesOrder;
      sortedOrder(originFiles, filesOrder);
      originFiles = fileNameList(originFiles, filesOrder);
      sortedOrder(restartFiles, filesOrder);
      restartFiles = fileNameList(restartFiles, filesOrder);
      // Check size of the files list
      if(!assertEquals(originFiles.size(), restartFiles.size(), &result)) {
        result.expected="originFiles.size() = " + name(originFiles.size());
        forAll(originFiles, i) {
          result.expected+="\n\t - " + originFiles[i];
        }
        result.expected+="\n";
        result.actual="restartFiles.size() = " + name(restartFiles.size());
        forAll(restartFiles, i) {
          result.actual+="\n\t - " + restartFiles[i];
        }
        result.actual+="\n";
        // Clean folder
        system(testFolder+"/Allclean > "+log_file+" 2>&1");
        system("rm -f "+testFolder+"/log.run");
        return result;
      }
      // Create time and mesh for origin and restard cases
      instant instantTime(originTime,originTimeDirs[dirI]);
      Foam::Time runTime_restart(restartTestFolder.path(),restartTestFolder.name());
      runTime_restart.setTime(instantTime,0);
      fvMesh mesh_restart
      (
          IOobject
          (
              regionName,
              runTime_restart.timeName(),
              runTime_restart,
              IOobject::MUST_READ
          )
      );
      Foam::Time runTime_origin(originTestFolder.path(),originTestFolder.name());
      runTime_origin.setTime(instantTime,0);
      fvMesh mesh_origin
      (
          IOobject
          (
              regionName,
              runTime_origin.timeName(),
              runTime_origin,
              IOobject::MUST_READ
          )
      );
      // For the PATOx boundary conditions
      simpleEnergyModel& energyModel_restart
      (
          meshLookupOrConstructModel<simpleEnergyModel>
          (
              mesh_restart,
              regionName,
              simpleEnergyModel::modelName
          )
      );
      simpleEnergyModel& energyModel_origin
      (
          meshLookupOrConstructModel<simpleEnergyModel>
          (
              mesh_origin,
              regionName,
              simpleEnergyModel::modelName
          )
      );
      // minimum value used in assertEquals function
      double min = 1e-14;
      // Check the files content
      forAll(restartFiles, fileI) {
        word fieldName = restartFiles[fileI].name();
        word className = findClass(restartFiles[fileI]);
        if (className=="volScalarField") {
          volScalarField* field_restart_ptr;
          volScalarField* field_origin_ptr;
          if (fieldName=="Ta") {
            field_restart_ptr=new volScalarField(energyModel_restart.refVolField<scalar>("Ta"));
            field_origin_ptr=new volScalarField(energyModel_origin.refVolField<scalar>("Ta"));
          } else {
            field_restart_ptr = new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh_restart.time().timeName(),
                    mesh_restart,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ), mesh_restart
            );
            field_origin_ptr = new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh_origin.time().timeName(),
                    mesh_origin,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ), mesh_restart
            );
          }
          volScalarField& field_restart = *field_restart_ptr;
          volScalarField& field_origin = *field_origin_ptr;
          if(!assertEquals(field_origin,field_restart,&result,min)) {
            // Clean folder
            system(testFolder+"/Allclean > "+log_file+" 2>&1");
            system("rm -f "+testFolder+"/log.run");
            return result;
          }
        } else {
          if (className=="volVectorField") {
            volVectorField field_restart
            (
                IOobject
                (
                    fieldName,
                    mesh_restart.time().timeName(),
                    mesh_restart,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ), mesh_restart
            );
            volVectorField field_origin
            (
                IOobject
                (
                    fieldName,
                    mesh_origin.time().timeName(),
                    mesh_origin,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ), mesh_restart
            );
            if(!assertEquals(field_origin,field_restart,&result,min)) {
              // Clean folder
              system(testFolder+"/Allclean > "+log_file+" 2>&1");
              system("rm -f "+testFolder+"/log.run");
              return result;
            }
          } else {
            if (className=="volTensorField") {
              volTensorField field_restart
              (
                  IOobject
                  (
                      fieldName,
                      mesh_restart.time().timeName(),
                      mesh_restart,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE
                  ), mesh_restart
              );
              volTensorField field_origin
              (
                  IOobject
                  (
                      fieldName,
                      mesh_origin.time().timeName(),
                      mesh_origin,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE
                  ), mesh_restart
              );
              if(!assertEquals(field_origin,field_restart,&result,min)) {
                // Clean folder
                system(testFolder+"/Allclean > "+log_file+" 2>&1");
                system("rm -f "+testFolder+"/log.run");
                return result;
              }
            } else {
              if (className=="volSymmTensorField") {
                volSymmTensorField field_restart
                (
                    IOobject
                    (
                        fieldName,
                        mesh_restart.time().timeName(),
                        mesh_restart,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ), mesh_restart
                );
                volSymmTensorField field_origin
                (
                    IOobject
                    (
                        fieldName,
                        mesh_origin.time().timeName(),
                        mesh_origin,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ), mesh_restart
                );
                if(!assertEquals(field_origin,field_restart,&result,min)) {
                  // Clean folder
                  system(testFolder+"/Allclean > "+log_file+" 2>&1");
                  system("rm -f "+testFolder+"/log.run");
                  return result;
                }
              } else {
                if (className=="surfaceScalarField" || className=="surfaceVectorField" || className=="surfaceTensorField"
                    || className=="pointScalarField" || className=="pointVectorField" || className=="pointTensorField"
                    || className=="scalarField" || className=="vectorField" || className=="tensorField" || className=="labelList") {
                  continue; // not implemented
                } else {
                  FatalError << "class name \"" << className << "\" not found" << exit(FatalError);
                }
              }
            }
          }
        }
      }
      index_restart++;
    }
    // Clean folder
    system(testFolder+"/Allclean > "+log_file+" 2>&1");
    system("rm -f "+testFolder+"/log.run");
    // Change the tolerance back for Linux
#ifndef __APPLE__
    SubTest::tol=old_tol;
#endif
    // Unmute Info
    Info.level=level_Info_;
    return result;
  }

  // Function to remove all spaces from a given string
  static string removeSpaces(string str) {
    List<int> indexes;
    for (int i = 0; str[i]; i++)
      if (str[i] != ' ')
        indexes.append(i);
    string new_str;
    for (int i = 0; i<indexes.size(); i++)
      new_str += str[indexes[i]];
    return new_str;
  }

  // Find the class in a OpenFOAM field file (e.g. volScalarField)
  static string findClass(fileName file) {
    fileName file_name = changeEnviVar(file);
    IFstream dataFile(file_name);
    if (!dataFile.good()) {
      FatalErrorInFunction
          << "Cannot read file " << file_name
          << exit(FatalError);
    }
    string class_str="class";
    while (dataFile) {
      string line;
      dataFile.getLine(line);
      int first = line.find(class_str);
      int last = line.find(";");
      if (first >= 0 && last >= 0) {
        string class_name = line.substr(first+class_str.size(),last-first-class_str.size());
        return removeSpaces(class_name);
      }
    }
    FatalError << "\"class\" not found in " << file << exit(FatalError);
    return "";
  }
};




