#include "subtest.h"
#include <vector>
#include "IOFunctions.H"

class TestIOFunctions : public SubTest
{
 public:

  TestIOFunctions() {
    testSuiteName = "TestIOFunctions";
    tests.push_back(test1);
    tests.push_back(test2);
    tests.push_back(testReadFileData);
    tests.push_back(testReadFileData_3D);
  }

  static TestResult test1() {
    std::string suiteName = "IO functions";
    std::string testName = "IO functions: readTecplotFileData";
    std::string testDescription = "Test readTecplotFileData POINT datapacking";
    TestResult result(suiteName, testName, 1, testDescription);
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    const fileName file = pato_dir + "/testsuites/global/IO/tecplot_point.dat";
    int info_level = Info.level;
    Info.level=0;
    const List<scalarList> data = readTecplotFileData(file);
    Info.level=info_level;
    const List<scalarList> expected_data = {{1,1,1},{2,3,4},{100,95,90},{50,50,50},{1,1,0.9}};
    if(!assertEquals(data.size(), expected_data.size() , &result)) {
      result.expected="data size expected = "+name(expected_data.size());
      result.actual="data size actual = "+name(data.size());
      return result;
    }
    forAll(expected_data, i) {
      if(!assertEquals(data[i].size(), expected_data[i].size() , &result)) {
        result.expected="data[i] size expected = "+name(expected_data[i].size());
        result.actual="data[i] size actual = "+name(data[i].size());
        return result;
      }
      forAll(expected_data[i], j) {
        scalar actualScalar_ = data[i][j];
        scalar expectedScalar_ = expected_data[i][j];
        if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
          result.expected=" expected = "+name(expectedScalar_);
          result.actual=" actual = "+name(actualScalar_);
          return result;
        }
      }
    }
    return result;
  }

  static TestResult test2() {
    std::string suiteName = "IO functions";
    std::string testName = "IO functions: readTecplotFileData";
    std::string testDescription = "Test readTecplotFileData BLOCK datapacking";
    TestResult result(suiteName, testName, 1, testDescription);
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    const fileName file = pato_dir + "/testsuites/global/IO/tecplot_block.dat";
    int info_level = Info.level;
    Info.level=0;
    const List<scalarList> data = readTecplotFileData(file);
    Info.level=info_level;
    const List<scalarList> expected_data = {{1,1,1},{2,3,4},{100,95,90},{50,50,50},{1,1,0.9}};
    if(!assertEquals(data.size(), expected_data.size() , &result)) {
      result.expected="data size expected = "+name(expected_data.size());
      result.actual="data size actual = "+name(data.size());
      return result;
    }
    forAll(expected_data, i) {
      if(!assertEquals(data[i].size(), expected_data[i].size() , &result)) {
        result.expected="data[i] size expected = "+name(expected_data[i].size());
        result.actual="data[i] size actual = "+name(data[i].size());
        return result;
      }
      forAll(expected_data[i], j) {
        scalar actualScalar_ = data[i][j];
        scalar expectedScalar_ = expected_data[i][j];
        if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
          result.expected=" expected = "+name(expectedScalar_);
          result.actual=" actual = "+name(actualScalar_);
          return result;
        }
      }
    }
    return result;
  }

  static TestResult testReadFileData() {
    std::string suiteName = "IO functions";
    std::string testName = "IO functions: readTecplotFileData";
    std::string testDescription = "Test readTecplotFileData";
    TestResult result(suiteName, testName, 1, testDescription);
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    const fileName file = pato_dir + "/testsuites/global/IO/file_data.dat";
    List<scalarList> expected_data= {{1,4},{2,5},{3,6}};
    List<scalarList> data=readFileData(file);

    forAll(expected_data, i) {
      if(!assertEquals(data[i].size(), expected_data[i].size() , &result)) {
        result.expected="data[i] size expected = "+name(expected_data[i].size());
        result.actual="data[i] size actual = "+name(data[i].size());
        return result;
      }
      forAll(expected_data[i], j) {
        scalar actualScalar_ = data[i][j];
        scalar expectedScalar_ = expected_data[i][j];
        if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
          result.expected=" expected = "+name(expectedScalar_);
          result.actual=" actual = "+name(actualScalar_);
          return result;
        }
      }
    }
    return result;
  }

  static TestResult testReadFileData_3D() {
    std::string suiteName = "IO functions";
    std::string testName = "IO functions: readTecplotFileData_3D";
    std::string testDescription = "Test readTecplotFileData 3D";
    TestResult result(suiteName, testName, 1, testDescription);
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    const fileName file = pato_dir + "/testsuites/global/IO/BoundaryConditions_ArcJet_cylinder_3D";
    List<scalarList> expected_data= {{0,0.1,40,40.1,80},{101325,101325,101325,101325,101325},{0.3e-2,0.3,0.3,0.3e-2,0.3e-2},{0,1e7,1e7,0,0},{1,1,1,1,0},{0.5,0.5,0.5,0.5,0.5},{300,300,300,300,300},{0,0,0,0,0}};
    List<scalarList> data=readFileData(file);

    forAll(expected_data, i) {
      if(!assertEquals(data[i].size(), expected_data[i].size() , &result)) {
        result.expected="data[i] size expected = "+name(expected_data[i].size());
        result.actual="data[i] size actual = "+name(data[i].size());
        return result;
      }
      forAll(expected_data[i], j) {
        scalar actualScalar_ = data[i][j];
        scalar expectedScalar_ = expected_data[i][j];
        if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
          result.expected=" expected = "+name(expectedScalar_);
          result.actual=" actual = "+name(actualScalar_);
          return result;
        }
      }
    }
    return result;
  }
};

