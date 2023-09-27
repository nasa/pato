
#include "PATOx.H"
#include "subtest.h"

#include <vector>


class TestChemistryModel : public SubTest
{
 public:

  TestChemistryModel() {
    testSuiteName = "TestChemistryModel";
    tests.push_back(test1);
  }

  static TestResult test1() {
    std::string suiteName = "TestChemistryModel";
    std::string testName = "TestChemistryModel: Test 1 - (couple word description)";
    std::string testDescription = "(Extended description if desired)";
    TestResult result(suiteName, testName, 1, testDescription);

    std::vector<double> m;

    if(!assertEquals((long)0,(long)m.size(), &result)) {
      return result;
    }

    return result;
  }

};
