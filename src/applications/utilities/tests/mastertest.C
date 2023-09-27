// INCLUDE NEW TEST SUITES HERE - DO NOT REMOVE THIS LINE
#include "testsuites/global/TestCodingStyle.C"
#include "testsuites/global/TestRestart.C"
#include "testsuites/MaterialModel/IOModel/TestIOModel.C"
#include "testsuites/MaterialModel/ChemistryModel/TestChemistryModel.C"
#include "testsuites/global/TestIOFunctions.C"
#include "testsuites/global/TestMathFunctions.C"
#include "testsuites/MaterialModel/BoundaryConditionsModel/TestBoundaryMapping.C"
#include "testsuites/MaterialModel/BoundaryConditionsModel/TestHeatFluxBC.C"
#include "testsuites/MaterialModel/MassModel/TestMassModel.C"
#include "testsuites/MaterialModel/MaterialPropertiesModel/TestMaterialPropertiesModel.C"
#include "testsuites/tutorials/TestTutorials_folders.C"
#include "testsuites/tutorials/TestTutorials.C"

#include "mastertest.h"

#include <iostream>

MasterTest::MasterTest()
{

//ADD NEW TEST SUITES HERE - DO NOT REMOVE THIS LINE
  testSuites.push_back(TestCodingStyle());
  testSuites.push_back(TestRestart());
  testSuites.push_back(TestIOModel());
  testSuites.push_back(TestChemistryModel());
  testSuites.push_back(TestIOFunctions());
  testSuites.push_back(TestMathFunctions());
  testSuites.push_back(TestBoundaryMapping());
  testSuites.push_back(TestHeatFluxBC());
  testSuites.push_back(TestMassModel());
  testSuites.push_back(TestMaterialPropertiesModel());
  testSuites.push_back(TestTutorials_folders());
  testSuites.push_back(TestTutorials());
}

std::vector<TestResult> MasterTest::runAllTests()
{

  std::vector<TestResult> failedTests;

  for( unsigned int i = 0; i < testSuites.size(); i++ ) {
    std::cout << "| Test Suite " << i << " : " << testSuites[i].testSuiteName << std::endl;
    std::vector<TestResult> testSuiteFailedTests = testSuites[i].runAllTests();
    failedTests.insert(failedTests.end(), testSuiteFailedTests.begin(), testSuiteFailedTests.end()) ;
  }

  return failedTests;

}

std::vector<TestResult> MasterTest::runTests(std::vector<int> indices)
{
  std::vector<TestResult> failedTests;

  for( unsigned int i = 0; i < indices.size(); i++ ) {
    std::cout << "| Test Suite " << indices[i] << " : " << testSuites[indices[i]].testSuiteName << std::endl;
    std::vector<TestResult> testSuiteFailedTests = testSuites[indices[i]].runAllTests();
    failedTests.insert(failedTests.end(), testSuiteFailedTests.begin(), testSuiteFailedTests.end()) ;
  }

  return failedTests;
}

void MasterTest::printTestSuites()
{
  std::cout << "To run all tests, run without arguments or with '-a'" << std::endl;
  std::cout << "To run subset of test suites, list suite indexes as \n    individual arguments or grouped (eg. '0-3')" << std::endl;
  std::cout << "Available Test Suites: " << std::endl;

  for(unsigned int i=0; i<testSuites.size(); i++) {
    std::cout << i << ": " << testSuites[i].testSuiteName << std::endl;
  }
}

int MasterTest::numTests()
{
  int total = 0;
  for( unsigned int i = 0; i < testSuites.size(); i++ ) {
    total += testSuites[i].numTests();
  }

  return total;
}

int MasterTest::numTests(std::vector<int> indices)
{
  int total = 0;
  for( unsigned int i = 0; i < indices.size(); i++ ) {
    total += testSuites[indices[i]].numTests();
  }

  return total;
}

bool MasterTest::validSubTests(std::vector<int> indices)
{
  for(unsigned int i=0; i<indices.size(); i++) {
    if(indices[i]<0||indices[i]>=(int)testSuites.size()) {
      return false;
    }
  }

  return true;
}
