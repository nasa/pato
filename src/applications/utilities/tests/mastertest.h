#ifndef MASTERTEST_H
#define MASTERTEST_H

#include <set>

#include "testframework/subtest.h"
#include "testframework/testresult.h"

class MasterTest
{
public:
    MasterTest();
    void printTestSuites();
    std::vector<TestResult> runAllTests();
    std::vector<TestResult> runTests(std::vector<int> indices);
    int numTests();
    int numTests(std::vector<int> indices);
    bool validSubTests(std::vector<int> indices);

private:
    std::vector<SubTest> testSuites;
};

#endif // MASTERTEST_H
