
#include "PATOx.H"
#include "subtest.h"

#include "fileOperation.H"
#include "IOFunctions.H"

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>


class TestIOModel : public SubTest
{
 public:

  TestIOModel() {
    testSuiteName = "TestIOModel";
    tests.push_back(test1);
//    tests.push_back(test2);
//        autoPtr<specifiedMassModel> newModel(specifiedMassModel::New(mesh));
  }

  static TestResult test1() {
    std::string suiteName = "TestIOModel";
    std::string testName = "TestIOModel: Test 1 - (couple word description)";
    std::string testDescription = "(Extended description if desired)";
    TestResult result(suiteName, testName, 1, testDescription);


    word caseName = "testWriteControl";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    //Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/IOModel",caseName);
    //word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/IOModel/" + caseName;

    // Run the case
    system(case_+"/Allclean > /dev/null");
    system(case_+"/Allrun > /dev/null");

    // Read all the directories
    fileNameList dirs=readDir(case_, fileType::directory);
    fileNameList tdirs;
    // Add all directories that are numbers and sort
    for (int i=0; i<dirs.size(); i++) {
      if (isDouble(dirs[i])) {
        tdirs.append(dirs[i]);
      }
    }
    std::sort(tdirs.begin(),tdirs.end());

    // Get written field names and sort (Capital letters are first)
    fileNameList vars=readDir(case_+"/"+tdirs[1]+"/porousMat",fileType::file);
    std::sort(vars.begin(), vars.end());

    // Get output times
    scalarList output = readFileData(case_+"/output/porousMat/scalar/Ta_plot")[0];

    // Clean the case
    system(case_+"/Allclean > /dev/null");

    // Check the number of directories and the times
    if(!assertEquals((long) 3,(long) tdirs.size(), &result)) {
      return result;
    }
    if(!assertEquals((double)0.203898, std::stod(tdirs[1]), &result)) {
      return result;
    }
    if(!assertEquals((double)0.701462, std::stod(tdirs[2]), &result)) {
      return result;
    }

    // Check that the files in the time directories are correct
    fileNameList correctVars;
    correctVars.append("Ta");
    correctVars.append("nu");
    correctVars.append("p");
    correctVars.append("rho_s");
    for (int i=0; i<vars.size(); i++) {
      if(!assertEquals((std::string)correctVars[i], (std::string)vars[i], &result)) {
        return result;
      }
    }

    // Check the output
    if(!assertEquals((long)4, (long)output.size(), &result)) {
      return result;
    }
    if(!assertEquals((double)0.402924, (double)output[1], &result)) {
      return result;
    }
    if(!assertEquals((double)0.601949, (double)output[2], &result)) {
      return result;
    }

    return result;
  }

  static TestResult test2() {
    std::string suiteName = "TestIOModel";
    std::string testName = "TestIOModel: Test 2 - (couple word description)";
    std::string testDescription = "(Extended description if desired)";
    TestResult result(suiteName, testName, 2, testDescription);

    std::vector<double> m;

    if(!assertEquals((long)0,(long)m.size(), &result)) {
      return result;
    }

    return result;
  }

  static bool isDouble(fileName myFile) {
    std::istringstream iss(myFile);
    double d;
    iss >> std::noskipws >> d; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail();
  }

};
