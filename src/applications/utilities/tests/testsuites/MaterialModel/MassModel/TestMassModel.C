
#include "PATOx.H"
#include "subtest.h"

#include <vector>


class TestMassModel : public SubTest
{
 public:

  TestMassModel() {
    testSuiteName = "TestMassModel";
    tests.push_back(test1);
    tests.push_back(test2);
//        autoPtr<specifiedMassModel> newModel(specifiedMassModel::New(mesh));
  }

  static TestResult test1() {
    std::string suiteName = "TestMassModel";
    std::string testName = "TestMassModel: Test 1 - (couple word description)";
    std::string testDescription = "(Extended description if desired)";
    TestResult result(suiteName, testName, 1, testDescription);

    std::vector<double> m;

    if(!assertEquals((long)0,(long)m.size(), &result)) {
      return result;
    }

    return result;
  }

  static TestResult test2() {
    std::string suiteName = "TestMassModel";
    std::string testName = "TestMassModel: Test 2 - (couple word description)";
    std::string testDescription = "(Extended description if desired)";
    TestResult result(suiteName, testName, 2, testDescription);

    std::vector<double> m;

    if(!assertEquals((long)0,(long)m.size(), &result)) {
      return result;
    }

    return result;
  }

};
