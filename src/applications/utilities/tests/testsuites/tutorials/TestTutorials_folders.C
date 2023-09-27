#include "subtest.h"
#include <vector>

class TestTutorials_folders : public SubTest
{
 public:

  TestTutorials_folders() {
    testSuiteName = "TestTutorials_folders";
    tests.push_back(checkTutorialFolders);
  }

  static TestResult checkTutorialFolders() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "Tutorials_folders";
    std::string testName = "Tutorials: Test 1 - Regression testing: inspect the tutorials and reference folders";
    std::string testDescription = "";
    TestResult result(suiteName, testName, 1, testDescription);
    fileName tutoFolder = getEnv("PATO_DIR")+"/tutorials";
    fileName refFolder = getEnv("PATO_DIR")+"/src/applications/utilities/tests/testsuites/tutorials/ref";
    fileNameList tutoFolders = searchFoldersKeyword(tutoFolder,"Allrun");
    fileNameList refFolders  = searchFoldersKeyword(refFolder,"Allrun");
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
    Info.level=level_Info_;
    return result;
  }
};
