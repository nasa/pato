#include "subtest.h"
#include <vector>
#include "mathFunctions.H"

class TestMathFunctions : public SubTest
{
 public:

  TestMathFunctions() {
    testSuiteName = "TestMathFunctions";
    tests.push_back(test1);
    tests.push_back(test2);
    tests.push_back(test3);
    tests.push_back(test4);
    tests.push_back(test5);
    tests.push_back(test6);
  }

  static TestResult test1() {
    std::string suiteName = "math functions";
    std::string testName = "math functions: bilinear";
    std::string testDescription = "Test bilinear interpolation";
    TestResult result(suiteName, testName, 1, testDescription);
    scalarList x;
    scalarList y;
    scalarList f_xy;
    x.append(0);
    x.append(0);
    x.append(0);
    x.append(1);
    x.append(1);
    x.append(1);
    y.append(1);
    y.append(2);
    y.append(3);
    y.append(3);
    y.append(2);
    y.append(1);
    f_xy.append(10);
    f_xy.append(20);
    f_xy.append(30);
    f_xy.append(30);
    f_xy.append(20);
    f_xy.append(10);
    scalar actualScalar_ = bilinearInterpolation(x,y,f_xy,0.5,2.5);
    scalar expectedScalar_=25;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" expected = "+name(expectedScalar_);
      result.actual=" actual = "+name(actualScalar_);
      return result;
    }
    return result;
  }

  static TestResult test2() {
    std::string suiteName = "math functions";
    std::string testName = "math functions: bilinear";
    std::string testDescription = "Test bilinear interpolation";
    TestResult result(suiteName, testName, 1, testDescription);
    scalarList x;
    scalarList y;
    scalarList f_xy;
    x.append(0);
    x.append(0);
    x.append(0);
    x.append(1);
    x.append(1);
    x.append(1);
    y.append(1);
    y.append(2);
    y.append(3);
    y.append(4);
    y.append(2);
    y.append(1);
    f_xy.append(10);
    f_xy.append(20);
    f_xy.append(30);
    f_xy.append(40);
    f_xy.append(20);
    f_xy.append(10);
    scalar actualScalar_ = bilinearInterpolation(x,y,f_xy,0.5,2.5);
    scalar expectedScalar_=25;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" expected = "+name(expectedScalar_);
      result.actual=" actual = "+name(actualScalar_);
      return result;
    }
    return result;
  }

  static TestResult test3() {
    std::string suiteName = "math functions";
    std::string testName = "math functions: bilinear";
    std::string testDescription = "Test bilinear interpolation";
    TestResult result(suiteName, testName, 1, testDescription);
    scalarList x;
    scalarList y;
    scalarList f_xy;
    x.append(0);
    x.append(0);
    x.append(0);
    x.append(1);
    x.append(1);
    x.append(1);
    y.append(1);
    y.append(2);
    y.append(3);
    y.append(4);
    y.append(2);
    y.append(1);
    f_xy.append(10);
    f_xy.append(20);
    f_xy.append(30);
    f_xy.append(40);
    f_xy.append(20);
    f_xy.append(10);
    scalar actualScalar_ = bilinearInterpolation(x,y,f_xy,1,4);
    scalar expectedScalar_=40;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" expected = "+name(expectedScalar_);
      result.actual=" actual = "+name(actualScalar_);
      return result;
    }
    return result;
  }

  static TestResult test4() {
    std::string suiteName = "math functions";
    std::string testName = "math functions: bilinear";
    std::string testDescription = "Test bilinear interpolation";
    TestResult result(suiteName, testName, 1, testDescription);
    scalarList x;
    scalarList y;
    scalarList f_xy;
    x.append(0);
    x.append(0);
    x.append(0);
    x.append(1);
    x.append(1);
    x.append(1);
    y.append(1);
    y.append(2);
    y.append(3);
    y.append(4);
    y.append(2);
    y.append(1);
    f_xy.append(10);
    f_xy.append(20);
    f_xy.append(30);
    f_xy.append(40);
    f_xy.append(20);
    f_xy.append(10);
    scalar actualScalar_ = bilinearInterpolation(x,y,f_xy,0,1);
    scalar expectedScalar_=10;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" expected = "+name(expectedScalar_);
      result.actual=" actual = "+name(actualScalar_);
      return result;
    }
    return result;
  }

  static TestResult test5() {
    std::string suiteName = "math functions";
    std::string testName = "math functions: bilinear";
    std::string testDescription = "Test bilinear interpolation";
    TestResult result(suiteName, testName, 1, testDescription);
    scalarList x;
    scalarList y;
    scalarList f_xy;
    x.append(0);
    x.append(0);
    x.append(0);
    x.append(1);
    x.append(1);
    x.append(1);
    y.append(1);
    y.append(2);
    y.append(3);
    y.append(4);
    y.append(2);
    y.append(1);
    f_xy.append(10);
    f_xy.append(20);
    f_xy.append(30);
    f_xy.append(40);
    f_xy.append(20);
    f_xy.append(10);
    scalar actualScalar_ = bilinearInterpolation(x,y,f_xy,0,3);
    scalar expectedScalar_=30;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" expected = "+name(expectedScalar_);
      result.actual=" actual = "+name(actualScalar_);
      return result;
    }
    return result;
  }

  static TestResult test6() {
    std::string suiteName = "math functions";
    std::string testName = "math functions: bilinear";
    std::string testDescription = "Test bilinear interpolation";
    TestResult result(suiteName, testName, 1, testDescription);
    scalarList x;
    scalarList y;
    scalarList f_xy;
    x.append(0);
    x.append(0);
    x.append(0);
    x.append(1);
    x.append(1);
    x.append(1);
    y.append(1);
    y.append(2);
    y.append(3);
    y.append(4);
    y.append(2);
    y.append(1);
    f_xy.append(10);
    f_xy.append(20);
    f_xy.append(30);
    f_xy.append(40);
    f_xy.append(20);
    f_xy.append(10);
    scalar actualScalar_ = bilinearInterpolation(x,y,f_xy,1,1);
    scalar expectedScalar_=10;
    if(!assertEquals(actualScalar_, expectedScalar_ , &result)) {
      result.expected=" expected = "+name(expectedScalar_);
      result.actual=" actual = "+name(actualScalar_);
      return result;
    }
    return result;
  }
};

