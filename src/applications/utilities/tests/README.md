# PATO Test

Testing Framework, developed for PATO software.

## Getting Started

The testing framework is based on qmake. To compile:  
1. Navigate to the build/ directory  
2. Execute "qmake" in a terminal to generate the make file  
3. Execute "make" to compile the code  

## Running Tests

To run all the testsuties, the executable can be run with no arguments, or -a. 

To run a subset of the testsuites, suite indices can be passed as arguments.   
   - Example: "./runtests 0 1 2"  
   - Example: "./runtests 0-2"  
   - Example: "./runtests 0-2 4"  

To see a list of available testsuites and their indices, run "./runtests -h" 

## Adding Custom Test Suite

To add a custom test suite:   
1. Create test suite class in a .cpp file in the testsuites/ directory   
    -use example_test.cpp   
2. Each test is defined by a static function, and a function pointer in the testsuite constructor   
3. Add the new testsuite class to the mastertest.cpp class  
    3a. Include at top of file (see example for example_test)  
    3b. Creat instance of testsuite class in the mastertest constructor (see example for example_test)  
4. Rerun qmake and make   

## Authors

* **Joseph Ferguson** - *joseph.c.ferguson@nasa.gov* 
* **Jeremie Meurisse** - *jeremie.b.meurisse@nasa.gov*



