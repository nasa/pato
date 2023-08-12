#include "subtest.h"
#include <vector>

class TestTutorials : public SubTest
{
 public:

  TestTutorials() {
    testSuiteName = "TestTutorials";
    tests.push_back(test_tuto_basic);
#if defined(FOAM_EXTEND)
    tests.push_back(test_tuto_foam_extend);
#endif
    tests.push_back(test_tuto_0D);
    tests.push_back(test_tuto_1D);
    tests.push_back(test_tuto_2D);
    tests.push_back(test_tuto_3D);
  }

  static TestResult test_tuto_foam_extend() {
    return test_tuto("foam-extend");
  }

  static TestResult test_tuto_basic() {
    return test_tuto("basic");
  }

  static TestResult test_tuto_0D() {
    return test_tuto("0D");
  }

  static TestResult test_tuto_1D() {
    return test_tuto("1D");
  }

  static TestResult test_tuto_2D() {
    return test_tuto("2D");
  }

  static TestResult test_tuto_3D() {
    return test_tuto("3D");
  }

  // Get the tutorial and reference folders and create the TestResult
  static TestResult test_tuto(word folder) {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string testName = "Tutorials: regression testing ("+folder+")";
    std::string testDescription = "Run the tutorials ("+folder+") and compare the results to the reference folder";
    const fileName pato_tuto_dir = getEnv("PATO_TUTORIALS");
    const fileName pato_unit_testing = getEnv("PATO_UNIT_TESTING");
    word reference_folder = "ref";
    createRefFolder_endTime_factor(reference_folder, folder); // Change the reference folder for endTime_factor != 1
    fileName tutoFolder = pato_tuto_dir+"/"+folder;
    fileName refFolder = pato_unit_testing+"/testsuites/tutorials/"+reference_folder+"/"+folder;
    fileNameList tutoFolders = searchFoldersKeyword(tutoFolder,"Allrun");
    fileNameList refFolders  = searchFoldersKeyword(refFolder,"Allrun");
    TestResult result(runTestTutorial(testName,testDescription,tutoFolders,refFolders));
    if (endTime_factor != 1) {
      word command = "rm -rf " + refFolder.path();
      system(command);
    }
    Info.level=level_Info_;
    return result;
  }

  // Change the reference folder for endTime_factor != 1
  static void createRefFolder_endTime_factor(word& reference_folder, word folder) {
    if (endTime_factor != 1) {
      const fileName pato_tuto_dir = getEnv("PATO_TUTORIALS");
      const fileName pato_unit_testing = getEnv("PATO_UNIT_TESTING");
      fileName tutoFolder = pato_tuto_dir+"/"+folder;
      fileName refFolder_origin = pato_unit_testing+"/testsuites/tutorials/"+reference_folder+"/"+folder;
      reference_folder+="_endTime_factor_"+name(endTime_factor);
      fileName refFolder = pato_unit_testing+"/testsuites/tutorials/"+reference_folder+"/"+folder;
      word command = "rm -rf " + refFolder + "; mkdir -p "+ refFolder.path() +"; cp -r " + refFolder_origin + " " + refFolder;
      system(command);
      fileNameList tutoFolders = searchFoldersKeyword(tutoFolder,"Allrun");
      fileNameList refFolders  = searchFoldersKeyword(refFolder,"Allrun");
      labelList visitOrder;
      sortedOrder(tutoFolders, visitOrder);
      tutoFolders = fileNameList(tutoFolders, visitOrder);
      sortedOrder(refFolders, visitOrder);
      refFolders = fileNameList(refFolders, visitOrder);

      if (tutoFolders.size() != refFolders.size()) {
        FatalErrorInFunction << "tutoFolders.size() != refFolders.size()" << exit(FatalError);
      }

      forAll(refFolders, i) {
        if (refFolders[i].name()!=tutoFolders[i].name()) {
          FatalErrorInFunction << refFolders[i] << " != " << tutoFolders[i] << exit(FatalError);
        }
        fileName controlDict_file = tutoFolders[i]+"/system/controlDict";
        if (!isFile(controlDict_file)) {
          FatalErrorInFunction << controlDict_file << " not found." << exit(FatalError);
        }
        word command = "echo $(grep \"endTime \\+[0-9]*;\" "+controlDict_file+" | tr -dc ‘0-9’)";
        scalar endTime = std::stof(exec(command.c_str()));
        command = "echo $(grep \"startTime \\+[0-9]*;\" "+controlDict_file+" | tr -dc ‘0-9’)";
        scalar startTime = std::stof(exec(command.c_str()));
        endTime=startTime+(endTime-startTime)/endTime_factor;

        fileNameList refFiles=filesInFolder(refFolders[i]+"/output");
        forAll(refFiles, j) {
          fileName file = refFiles[j];
          if (file.name()=="empty") {
            continue;
          }
          if (file.path().name()=="profile") {
            string file_name = file.name();
            size_t last_index = file_name.find_last_not_of("0123456789");
            word time_name = file_name.substr(last_index + 1);
            if (isNumber(time_name)) {
              scalar time = stof(time_name);
              if (time > endTime) {
                fileName file_copy = file;
                file_copy.replaceAll("(","\\(");
                file_copy.replaceAll(")","\\)");
                command = "rm -f " + file_copy;
                system(command);
              }
            } else {
              changeFile_endTime_factor(file, endTime);
            }
          } else {
            changeFile_endTime_factor(file, endTime);
          }
        }
      }
    }
  }

  // Change the files in the reference folders for endTime_factor != 1
  static void changeFile_endTime_factor(fileName file, scalar endTime) {
    List<scalarList> data = readFileData(file);
    int index = -1;
    forAll(data[0], dataJ) {
      if (data[0][dataJ]>endTime) {
        index = dataJ;
        if (dataJ == 0) {
          FatalErrorInFunction << file << ": data[0][0]("<<data[0][0]<<")>endTime(" << endTime << ")" << exit(FatalError);
        }
        break;
      }
    }
    fileName file_copy = file;
    file_copy.replaceAll("(","\\(");
    file_copy.replaceAll(")","\\)");
    word command = "rm -f " + file_copy;
    OFstream os(file);
    forAll(data[0], dataJ) {
      if (dataJ >= index) {
        break;
      }
      forAll(data, dataI) {
        os << data[dataI][dataJ] << " ";
      }
      os << endl;
    }
  }

  // Modify the controlDict in the tutorial folder
  static void modifyControlDict(fileName tutoFolder) {
    if (endTime_factor != 1) {
      fileName controlDict_file = tutoFolder+"/system/controlDict";
      if (!isFile(controlDict_file)) {
        FatalErrorInFunction << controlDict_file << " not found." << exit(FatalError);
      }
      word command = "echo $(grep \"endTime \\+[0-9]*;\" "+controlDict_file+" | tr -dc ‘0-9’)";
      scalar endTime = std::stof(exec(command.c_str()));
      command = "echo $(grep \"startTime \\+[0-9]*;\" "+controlDict_file+" | tr -dc ‘0-9’)";
      scalar startTime = std::stof(exec(command.c_str()));
      endTime=startTime+(endTime-startTime)/endTime_factor;
      command = "cp " + controlDict_file + " " + controlDict_file + ".copy";
      system(command);
      command = sed_command + " -i \"s/endTime \\+[0-9]*;/endTime "+ name(endTime) +";/g\" " + controlDict_file;
      system(command);
    }

  }

  // Put back the controlDict in the tutorial folder
  static void deleteControlDictCopy(fileName tutoFolder) {
    if (endTime_factor != 1) {
      fileName controlDict_file = tutoFolder+"/system/controlDict";
      if (!isFile(controlDict_file+".copy")) {
        FatalErrorInFunction << controlDict_file << ".copy not found." << exit(FatalError);
      }
      word command = "cp " + controlDict_file + ".copy " + controlDict_file + "; rm -f " + controlDict_file + ".copy ";
      system(command);
    }
  }

  // Run the tutorial folders and compare them to the reference folders
  static TestResult runTestTutorial(std::string testName, std::string testDescription, fileNameList tutoFolders, fileNameList refFolders) {
    std::string suiteName = "Tutorials";
    TestResult result(suiteName, testName, 1, testDescription);

    labelList visitOrder;
    sortedOrder(tutoFolders, visitOrder);
    tutoFolders = fileNameList(tutoFolders, visitOrder);
    sortedOrder(refFolders, visitOrder);
    refFolders = fileNameList(refFolders, visitOrder);

    if(!assertEquals(refFolders.size(),tutoFolders.size(), &result)) {
      result.expected="refFolders.size() = " + name(refFolders.size());
      forAll(refFolders, i) {
        result.expected+="\n\t - " + refFolders[i];
      }
      result.expected+="\n";
      result.actual="tutoFolders.size() = " + name(tutoFolders.size());
      forAll(tutoFolders, i) {
        result.actual+="\n\t - " + tutoFolders[i];
      }
      result.actual+="\n";
      return result;
    }

    forAll(tutoFolders,folderI) {
      fileName folder = tutoFolders[folderI];
      fileName reffolder = refFolders[folderI];
      const fileName pato_tuto =  getEnv("PATO_TUTORIALS");
      Info.level = 2;
      Info << "\t - [" << folderI+1 << "/" << tutoFolders.size() << "] " << (word) testName << " : " << elapsedClockTime() << " - " << tutoFolders[folderI] << endl;
      Info.level = 0;
      fileName log_file = "log.runTestTutorial.";
      fileName log_folder = folder;
      log_folder.replaceAll(pato_tuto,"");
      forAll(log_folder.components(), i) {
        log_file+=log_folder.components()[i];
        if ( i < log_folder.components().size()-1) {
          log_file+="_";
        }
      }

      modifyControlDict(folder); // modify the endTime in the controlDict files
      system(folder+"/Allclean > "+log_file);
      system(folder+"/Allrun > "+log_file);
      deleteControlDictCopy(folder); // remove the controlDict.copy files

      fileNameList tutoFiles=filesInFolder(folder+"/output");
      fileNameList refFiles=filesInFolder(reffolder+"/output");

      labelList visitOrder;
      sortedOrder(tutoFiles, visitOrder);
      tutoFiles = fileNameList(tutoFiles, visitOrder);
      sortedOrder(refFiles, visitOrder);
      refFiles = fileNameList(refFiles, visitOrder);

      if(!assertEquals(refFiles.size(), tutoFiles.size(), &result)) {
        result.expected="refFiles.size() = " + name(refFiles.size());
        forAll(refFiles, i) {
          result.expected+="\n\t - " + refFiles[i];
        }
        result.expected+="\n";
        result.actual="tutoFiles.size() = " + name(tutoFiles.size());
        forAll(tutoFiles, i) {
          result.actual+="\n\t - " + tutoFiles[i];
        }
        result.actual+="\n";
        return result;
      }

      forAll(tutoFiles, fileI) {
        fileName tutofile = tutoFiles[fileI];
        fileName reffile = refFiles[fileI];
        if(!assertEquals(tutofile,reffile,&result)) {
          return result;
        }
      }
      system(folder+"/Allclean > "+log_file);
    }
    system("rm log.runTestTutorial.*");
    return result;
  }

};
