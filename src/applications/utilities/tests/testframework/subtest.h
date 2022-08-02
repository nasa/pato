#ifndef SUBTEST_H
#define SUBTEST_H

#include <vector>
#include <string>
#include <cmath>

#include "testresult.h"
#include "fvCFD.H"
#include "IOFunctions.H"

class SubTest;
typedef TestResult (*test)();

class SubTest
{

public:

    SubTest();
    std::vector<TestResult> runAllTests();
    int numTests();
    void addTest(test t);
    std::string testSuiteName;

protected:

    std::vector<test> tests;
    static List<scalar> tol; // tolerances for double/float
    static scalar endTime_factor; // endTime factor: endTime = startTime + (endTime - startTime) / endTime_factor;
    static word sed_command; // sed command
    static Foam::clock clock_; // clock object with the start time
    static word elapsedClockTime(); // print the elapsed clock time
    static word exec(const char* cmd); // execute a command and gets the results
    static void init_options();  // initialize the endTime factor and the tolerances for double/float
    static bool assertEquals(int valExpected, int valActual, TestResult *result);
    static bool assertEquals(long valExpected, long valActual, TestResult *result);
    static bool assertEquals(float valExpected, float valActual, TestResult *result);
    static bool assertEquals(double valExpected, double valActual, TestResult *result);
    static bool assertEquals(std::string valExpected, std::string valActual, TestResult *result);
    static bool assertEquals(char valExpected, char valActual, TestResult *result);
    static bool assertEquals(bool valExpected, bool valActual, TestResult *result);
    static bool assertEquals(vector valExpected, vector valActual, TestResult *result);
    static bool assertEquals(tensor valExpected, tensor valActual, TestResult *result);
    static bool assertEquals(fileName valExpected, fileName valActual, TestResult *result);
};

#endif // SUBTEST_H
