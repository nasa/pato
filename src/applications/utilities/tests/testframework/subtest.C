#include "subtest.h"

#include <iostream>

#include "toString.h"
#include "color.h"
#include <iomanip>      // std::setprecision
#include <sstream> // stringstream



List<scalar> SubTest::tol;
scalar SubTest::endTime_factor;
word SubTest::sed_command;
Foam::clock SubTest::clock_;

SubTest::SubTest()
{
  testSuiteName = "";
  sed_command = "sed";
#ifdef __APPLE__
  sed_command = "gsed";
#endif
  init_options();
}

Foam::word SubTest::elapsedClockTime()
{
  std::ostringstream osBuffer;
  time_t t = clock_.elapsedClockTime();
  struct tm *timeStruct = gmtime(&t);
  osBuffer
      << std::setfill('0') << "ClockTime = "
      << std::setw(2) << timeStruct->tm_hour
      << 'h' << std::setw(2) << timeStruct->tm_min
      << 'm' << std::setw(2) << timeStruct->tm_sec
      << 's';
  return osBuffer.str();
}

Foam::word SubTest::exec(const char* cmd)
{
  char buffer[128];
  word result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (fgets(buffer, sizeof buffer, pipe) != NULL) {
      result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return result;
}

void SubTest::init_options()
{
  int n_tol = 2;
#ifdef __APPLE__
  n_tol = 4;
#endif
  if (tol.size() != n_tol) {
    tol.resize(n_tol);
    const fileName pato_test = getEnv("PATO_UNIT_TESTING");
    fileName file = pato_test+"/testframework/runtests_options";
    Foam::Time runTime(file.path(),file.name());
    word handlerType = fileOperation::defaultFileHandler;
    int level = Info.level;
    Info.level=0;
    autoPtr<fileOperation> handler
    (
        fileOperation::New
        (
            handlerType,
            writeInfoHeader
        )
    );
    Info.level=level;
    Foam::fileHandler(handler);
    IOdictionary dict(
        IOobject
        (
            "",
            "",
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    forAll(tol, i) {
      scalar value = 0;
#ifdef __APPLE__
      value = readScalar(dict.lookup("tol_mac_"+name(i)));
#else
      value = readScalar(dict.lookup("tol_linux_"+name(i)));
#endif
      tol[i]=value;
    }
    endTime_factor = readScalar(dict.lookup("endTime_factor"));
  }
}

int SubTest::numTests()
{
  return tests.size();
}

void SubTest::addTest(test t)
{
  tests.push_back(t);
}

std::vector<TestResult> SubTest::runAllTests()
{

  Color::Modifier red(Color::FG_RED);
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);

  std::vector<TestResult> failedTests;

  for(unsigned int i=0; i<tests.size(); i++) {

    TestResult result = (tests[i])();

    if(result.failed == 0) {
      std::cout << "[-- " << testSuiteName << " -- Test " << i+1 << " / " << tests.size() <<  ": " << green << "Success" << def << " -- " << elapsedClockTime() << " --]"  << std::endl;
    } else {
      failedTests.push_back(result);
      std::cout << red << "[-- " << testSuiteName << " -- Test " << i+1 << " / " << tests.size() <<  ": FAILED " << " -- " << elapsedClockTime() << " --]" << std::endl;
      std::cout << "   [-- Expected: " << result.expected << " -- Actual: " << result.actual << " --]" << def << std::endl;
    }
  }

  return failedTests;
}

bool SubTest::assertEquals(const volScalarField& expectedField, const volScalarField& actualField, TestResult *result)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    if(!assertEquals(expectedField[cellI],actualField[cellI],result)) {
      result->failed=1;
      result->expected=expectedField.name()+".bf()["+name(cellI)+"]="+name(expectedField[cellI]);
      result->actual=actualField.name()+".bf()["+name(cellI)+"]="+name(actualField[cellI]);
      return !result->failed; // failed
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      if(!assertEquals(expectedField.boundaryField()[patchI][faceI],actualField.boundaryField()[patchI][faceI],result)) {
        result->failed=1;
        result->expected=expectedField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI]);
        result->actual=actualField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI]);
        return !result->failed; // failed
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volVectorField& expectedField, const volVectorField& actualField, TestResult *result)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    for (int i = 0; i < 3 ; i++) {
      if(!assertEquals(expectedField[cellI][i],actualField[cellI][i],result)) {
        result->failed=1;
        result->expected=expectedField.name()+"["+name(cellI)+"]="+name(expectedField[cellI][i]);
        result->actual=actualField.name()+"["+name(cellI)+"]="+name(actualField[cellI][i]);
        return !result->failed; // failed
      }
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      for (int i = 0; i < 3 ; i++) {
        if(!assertEquals(expectedField.boundaryField()[patchI][faceI][i],actualField.boundaryField()[patchI][faceI][i],result)) {
          result->failed=1;
          result->expected=expectedField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI][i]);
          result->actual=actualField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI][i]);
          return !result->failed; // failed
        }
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volSymmTensorField& expectedField, const volSymmTensorField& actualField, TestResult *result)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    for (int i = 0; i < 6 ; i++) {
      if(!assertEquals(expectedField[cellI][i],actualField[cellI][i],result)) {
        result->failed=1;
        result->expected=expectedField.name()+"["+name(cellI)+"]="+name(expectedField[cellI][i]);
        result->actual=actualField.name()+"["+name(cellI)+"]="+name(actualField[cellI][i]);
        return !result->failed; // failed
      }
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      for (int i = 0; i < 6 ; i++) {
        if(!assertEquals(expectedField.boundaryField()[patchI][faceI][i],actualField.boundaryField()[patchI][faceI][i],result)) {
          result->failed=1;
          result->expected=expectedField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI][i]);
          result->actual=actualField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI][i]);
          return !result->failed; // failed
        }
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volTensorField& expectedField, const volTensorField& actualField, TestResult *result)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    for (int i = 0; i < 9 ; i++) {
      if(!assertEquals(expectedField[cellI][i],actualField[cellI][i],result)) {
        result->failed=1;
        result->expected=expectedField.name()+"["+name(cellI)+"]="+name(expectedField[cellI][i]);
        result->actual=actualField.name()+"["+name(cellI)+"]="+name(actualField[cellI][i]);
        return !result->failed; // failed
      }
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      for (int i = 0; i < 9 ; i++) {
        if(!assertEquals(expectedField.boundaryField()[patchI][faceI][i],actualField.boundaryField()[patchI][faceI][i],result)) {
          result->failed=1;
          result->expected=expectedField.name()+"["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI][i]);
          result->actual=actualField.name()+"["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI][i]);
          return !result->failed; // failed
        }
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volScalarField& expectedField, const volScalarField& actualField, TestResult *result, double min)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    if(!assertEquals(expectedField[cellI],actualField[cellI],result,min)) {
      result->failed=1;
      result->expected=expectedField.name()+".bf()["+name(cellI)+"]="+name(expectedField[cellI]);
      result->actual=actualField.name()+".bf()["+name(cellI)+"]="+name(actualField[cellI]);
      return !result->failed; // failed
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      if(!assertEquals(expectedField.boundaryField()[patchI][faceI],actualField.boundaryField()[patchI][faceI],result,min)) {
        result->failed=1;
        result->expected=expectedField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI]);
        result->actual=actualField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI]);
        return !result->failed; // failed
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volVectorField& expectedField, const volVectorField& actualField, TestResult *result, double min)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    for (int i = 0; i < 3 ; i++) {
      if(!assertEquals(expectedField[cellI][i],actualField[cellI][i],result,min)) {
        result->failed=1;
        result->expected=expectedField.name()+"["+name(cellI)+"]="+name(expectedField[cellI][i]);
        result->actual=actualField.name()+"["+name(cellI)+"]="+name(actualField[cellI][i]);
        return !result->failed; // failed
      }
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      for (int i = 0; i < 3 ; i++) {
        if(!assertEquals(expectedField.boundaryField()[patchI][faceI][i],actualField.boundaryField()[patchI][faceI][i],result,min)) {
          result->failed=1;
          result->expected=expectedField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI][i]);
          result->actual=actualField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI][i]);
          return !result->failed; // failed
        }
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volSymmTensorField& expectedField, const volSymmTensorField& actualField, TestResult *result, double min)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    for (int i = 0; i < 6 ; i++) {
      if(!assertEquals(expectedField[cellI][i],actualField[cellI][i],result,min)) {
        result->failed=1;
        result->expected=expectedField.name()+"["+name(cellI)+"]="+name(expectedField[cellI][i]);
        result->actual=actualField.name()+"["+name(cellI)+"]="+name(actualField[cellI][i]);
        return !result->failed; // failed
      }
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      for (int i = 0; i < 6 ; i++) {
        if(!assertEquals(expectedField.boundaryField()[patchI][faceI][i],actualField.boundaryField()[patchI][faceI][i],result,min)) {
          result->failed=1;
          result->expected=expectedField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI][i]);
          result->actual=actualField.name()+".bf()["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI][i]);
          return !result->failed; // failed
        }
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(const volTensorField& expectedField, const volTensorField& actualField, TestResult *result, double min)
{
  result->failed=0; // successed
  // Compare the cell-center values
  forAll(expectedField, cellI) {
    for (int i = 0; i < 9 ; i++) {
      if(!assertEquals(expectedField[cellI][i],actualField[cellI][i],result,min)) {
        result->failed=1;
        result->expected=expectedField.name()+"["+name(cellI)+"]="+name(expectedField[cellI][i]);
        result->actual=actualField.name()+"["+name(cellI)+"]="+name(actualField[cellI][i]);
        return !result->failed; // failed
      }
    }
  }
  // Compare the boundary values
  forAll(expectedField.boundaryField(), patchI) {
    forAll(expectedField.boundaryField()[patchI], faceI) {
      for (int i = 0; i < 9 ; i++) {
        if(!assertEquals(expectedField.boundaryField()[patchI][faceI][i],actualField.boundaryField()[patchI][faceI][i],result,min)) {
          result->failed=1;
          result->expected=expectedField.name()+"["+name(patchI)+"]["+name(faceI)+"]="+name(expectedField.boundaryField()[patchI][faceI][i]);
          result->actual=actualField.name()+"["+name(patchI)+"]["+name(faceI)+"]="+name(actualField.boundaryField()[patchI][faceI][i]);
          return !result->failed; // failed
        }
      }
    }
  }
  return !result->failed; // successed
}

bool SubTest::assertEquals(int valExpected, int valActual, TestResult *result)
{

  if(valExpected==valActual) {
    result->failed=0;
  } else {
    result->failed=1;
  }

  result->expected=toString::convert(valExpected);
  result->actual=toString::convert(valActual);

  return !result->failed;
}

bool SubTest::assertEquals(long valExpected, long valActual, TestResult *result)
{
  // If val < machine precision, val = 0
  double min = 1e-15;
  if(fabs(valActual) <= min) {
    valActual=0;
  }
  if(valActual==0||valExpected==0) {
    if(fabs(valExpected-valActual)<tol[0]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
  } else {
#ifdef __APPLE__
    if(valExpected < tol[1] || valActual < tol[1]) {
      if(fabs(1-valExpected/valActual)<tol[2]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    } else {
      if(fabs(1-valExpected/valActual)<tol[3]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    }
#else
    if(fabs(1-valExpected/valActual)<tol[1]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
#endif
  }

  std::stringstream stream_expected;
  stream_expected << std::fixed << std::setprecision(15) << valExpected;
  std::string s_expected = stream_expected.str();
  result->expected= s_expected;

  std::stringstream stream_actual;
  stream_actual << std::fixed << std::setprecision(15) << valActual;
  std::string s_actual = stream_actual.str();
  result->actual= s_actual;

  return !result->failed;
}

bool SubTest::assertEquals(float valExpected, float valActual, TestResult *result)
{
  // If val < machine precision, val = 0
  double min = 1e-15;
  if(fabs(valActual) <= min) {
    valActual=0;
  }
  if(valActual==0||valExpected==0) {
    if(fabs(valExpected-valActual)<tol[0]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
  } else {
#ifdef __APPLE__
    if(valExpected < tol[1] || valActual < tol[1]) {
      if(fabs(1-valExpected/valActual)<tol[2]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    } else {
      if(fabs(1-valExpected/valActual)<tol[3]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    }
#else
    if(fabs(1-valExpected/valActual)<tol[1]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
#endif
  }

  std::stringstream stream_expected;
  stream_expected << std::fixed << std::setprecision(15) << valExpected;
  std::string s_expected = stream_expected.str();
  result->expected= s_expected;

  std::stringstream stream_actual;
  stream_actual << std::fixed << std::setprecision(15) << valActual;
  std::string s_actual = stream_actual.str();
  result->actual= s_actual;

  return !result->failed;
}

bool SubTest::assertEquals(double valExpected, double valActual, TestResult *result)
{
  double min = 1e-15;
  if(fabs(valActual) <= min) {
    valActual=0;
  }
  if(valActual==0||valExpected==0) {
    if(fabs(valExpected-valActual)<tol[0]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
  } else {
#ifdef __APPLE__
    if(valExpected < tol[1] || valActual < tol[1]) {
      if(fabs(1-valExpected/valActual)<tol[2]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    } else {
      if(fabs(1-valExpected/valActual)<tol[3]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    }
#else
    if(fabs(1-valExpected/valActual)<tol[1]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
#endif
  }

  std::stringstream stream_expected;
  stream_expected << std::fixed << std::setprecision(15) << valExpected;
  std::string s_expected = stream_expected.str();
  result->expected= s_expected;

  std::stringstream stream_actual;
  stream_actual << std::fixed << std::setprecision(15) << valActual;
  std::string s_actual = stream_actual.str();

  result->actual= s_actual;

  return !result->failed;
}

bool SubTest::assertEquals(long valExpected, long valActual, TestResult *result, double min)
{
  // If val < machine precision, val = 0
  if(fabs(valActual) <= min) {
    valActual=0;
  }
  if(valActual==0||valExpected==0) {
    if(fabs(valExpected-valActual)<tol[0]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
  } else {
#ifdef __APPLE__
    if(valExpected < tol[1] || valActual < tol[1]) {
      if(fabs(1-valExpected/valActual)<tol[2]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    } else {
      if(fabs(1-valExpected/valActual)<tol[3]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    }
#else
    if(fabs(1-valExpected/valActual)<tol[1]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
#endif
  }

  std::stringstream stream_expected;
  stream_expected << std::fixed << std::setprecision(15) << valExpected;
  std::string s_expected = stream_expected.str();
  result->expected= s_expected;

  std::stringstream stream_actual;
  stream_actual << std::fixed << std::setprecision(15) << valActual;
  std::string s_actual = stream_actual.str();
  result->actual= s_actual;

  return !result->failed;
}

bool SubTest::assertEquals(float valExpected, float valActual, TestResult *result, double min)
{
  // If val < machine precision, val = 0
  if(fabs(valActual) <= min) {
    valActual=0;
  }
  if(valActual==0||valExpected==0) {
    if(fabs(valExpected-valActual)<tol[0]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
  } else {
#ifdef __APPLE__
    if(valExpected < tol[1] || valActual < tol[1]) {
      if(fabs(1-valExpected/valActual)<tol[2]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    } else {
      if(fabs(1-valExpected/valActual)<tol[3]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    }
#else
    if(fabs(1-valExpected/valActual)<tol[1]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
#endif
  }

  std::stringstream stream_expected;
  stream_expected << std::fixed << std::setprecision(15) << valExpected;
  std::string s_expected = stream_expected.str();
  result->expected= s_expected;

  std::stringstream stream_actual;
  stream_actual << std::fixed << std::setprecision(15) << valActual;
  std::string s_actual = stream_actual.str();
  result->actual= s_actual;

  return !result->failed;
}

bool SubTest::assertEquals(double valExpected, double valActual, TestResult *result, double min)
{
  if(fabs(valActual) <= min) {
    valActual=0;
  }
  if(valActual==0||valExpected==0) {
    if(fabs(valExpected-valActual)<tol[0]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
  } else {
#ifdef __APPLE__
    if(valExpected < tol[1] || valActual < tol[1]) {
      if(fabs(1-valExpected/valActual)<tol[2]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    } else {
      if(fabs(1-valExpected/valActual)<tol[3]) {
        result->failed=0;
      } else {
        result->failed=1;
      }
    }
#else
    if(fabs(1-valExpected/valActual)<tol[1]) {
      result->failed=0;
    } else {
      result->failed=1;
    }
#endif
  }

  std::stringstream stream_expected;
  stream_expected << std::fixed << std::setprecision(15) << valExpected;
  std::string s_expected = stream_expected.str();
  result->expected= s_expected;

  std::stringstream stream_actual;
  stream_actual << std::fixed << std::setprecision(15) << valActual;
  std::string s_actual = stream_actual.str();

  result->actual= s_actual;

  return !result->failed;
}

bool SubTest::assertEquals(std::string valExpected, std::string valActual, TestResult *result)
{

  if(valExpected==valActual) {
    result->failed=0;
  } else {
    result->failed=1;
  }
  result->expected=valExpected;
  result->actual=valActual;

  return !result->failed;
}

bool SubTest::assertEquals(char valExpected, char valActual, TestResult *result)
{

  if(valExpected==valActual) {
    result->failed=0;
  } else {
    result->failed=1;
  }
  result->expected=toString::convert(valExpected);
  result->actual=toString::convert(valActual);

  return !result->failed;
}

bool SubTest::assertEquals(bool valExpected, bool valActual, TestResult *result)
{

  if(valExpected==valActual) {
    result->failed=0;
  } else {
    result->failed=1;
  }
  result->expected=toString::convert(valExpected);
  result->actual=toString::convert(valActual);

  return !result->failed;
}

bool SubTest::assertEquals(vector valExpected, vector valActual, TestResult *result)
{

  if(valExpected==valActual) {
    result->failed=0;
  } else {
    result->failed=1;
  }
  result->expected="vector = " + name(valExpected) ;
  result->actual="vector = " + name(valActual);

  return !result->failed;
}

bool SubTest::assertEquals(tensor valExpected, tensor valActual, TestResult *result)
{

  if(valExpected==valActual) {
    result->failed=0;
  } else {
    result->failed=1;
  }
  result->expected="tensor = " + name(valExpected) ;
  result->actual="tensor = " + name(valActual);

  return !result->failed;
}


bool SubTest::assertEquals(fileName valExpected, fileName valActual, TestResult *result)
{
  List<scalarList> data_expected(readFileData(valExpected));
  List<scalarList> data_actual(readFileData(valActual));

  if(!assertEquals(data_expected.size(),data_actual.size(),result)) {
    result->expected="\nfileName = " + valExpected + "\n column size = " + name(data_expected.size()) + "\n";
    result->actual="\nfileName = " + valActual + "\n column size = " + name(data_actual.size()) + "\n";
    return !result->failed;
  }

  forAll(data_expected, dataI) {
    if(!assertEquals(data_expected[dataI].size(),data_actual[dataI].size(),result)) {
      result->expected="\nfileName = " + valExpected + "\n row size = " + name(data_expected[dataI].size()) + "\n";
      result->actual="\nfileName = " + valActual + "\n row size = " + name(data_actual[dataI].size())+ "\n";
      return !result->failed;
    }
    forAll(data_expected[dataI],dataJ) {
      if(!assertEquals(data_expected[dataI][dataJ],data_actual[dataI][dataJ],result)) {
        result->expected="\nfileName = " + valExpected + "\n data["+name(dataJ)+","+name(dataI)+"] = " + name(data_expected[dataI][dataJ]) + "\n";
        result->actual="\nfileName = " + valActual + "\n data["+name(dataJ)+","+name(dataI)+"] = " + name(data_actual[dataI][dataJ])+ "\n";
        return !result->failed;
      }
    }
  }

  return !result->failed;
}
