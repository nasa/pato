
#include "PATOx.H"
#include "subtest.h"

#include "fileOperation.H"
#include "IOFunctions.H"

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>


class TestHeatFluxBC : public SubTest
{
 public:

  TestHeatFluxBC() {
    testSuiteName = "TestHeatFluxBC";
    tests.push_back(test1);
  }

  static TestResult test1() {
    std::string suiteName = "TestHeatFluxBC";
    std::string testName = "TestHeatFluxBC: Test 1 - MSL Test";
    std::string testDescription = "(Extended description if desired)";
    TestResult result(suiteName, testName, 1, testDescription);


    word caseName = "testHeatFluxBC";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    //Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/IOModel",caseName);
    //word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    // Run the case
    system(case_+"/Allclean > /dev/null");
    system(case_+"/Allrun > /dev/null");

    // Get output times
    scalarList output = readFileData(case_+"/output/porousMat/scalar/Ta_surfacePatch")[1];

    // Clean the case
    system(case_+"/Allclean > /dev/null");

    // Check the output
    if(!assertEquals((long)2, (long)output.size(), &result)) {
      return result;
    }
    if(!assertEquals((double)719.095055184, (double)output[1], &result)) {
      return result;
    }

    return result;
  }

};
