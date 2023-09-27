#include "subtest.h"
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

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

class TestCodingStyle : public SubTest
{
 public:

  TestCodingStyle() {
    testSuiteName = "TestCodingStyle";
    tests.push_back(test);
  }

  static TestResult test() {
    std::string suiteName = "PATO Coding Style";
    std::string testName = "PATO Coding Style in $PATO_DIR";
    std::string testDescription = "Test the coding style in $PATO_DIR using patoCodingStyle utility";
    TestResult result(suiteName, testName, 1, testDescription);
    word testing_dir = "$PATO_DIR";
    word cmd = ". $PATO_DIR/src/applications/utilities/runFunctions/RunFunctions; pato_init; patoCodingStyle -dir "
               +testing_dir;
    int level_Info = Info.level;
    Info.level = 0;
    string output = exec(cmd.c_str());
    Info.level = level_Info;
    int found_no_form = output.find("No formatted files");
    int actualScalar_ = -1;
    if (found_no_form >= 0) {
      actualScalar_=0;
    } else {
      int found_form = output.find("Total number of formatted files = ");
      if (found_form >= 0) {
        int size_str = string("Total number of formatted files = ").size();
        string num_form_files = output.substr(found_form+size_str);
        actualScalar_=stoi(num_form_files);
      }
    }
    int expectedScalar_= 0;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" Number of formatted files in $PATO_DIR = "+name(expectedScalar_);
      result.actual=" Number of formatted files in $PATO_DIR = "+name(actualScalar_)+
                    " --\n                  Please run \"patoCodingStyle -dir $PATO_DIR\" before to commit";
      return result;
    }
    return result;
  }
};

