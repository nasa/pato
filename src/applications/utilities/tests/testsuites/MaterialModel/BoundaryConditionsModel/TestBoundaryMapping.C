#include "subtest.h"
#include <vector>
#include "simpleBoundaryMappingModel.H"
#include <iomanip>      // std::setprecision
#include "IOmanip.H"

class TestBoundaryMapping : public SubTest
{
 public:

  TestBoundaryMapping() {
    testSuiteName = "TestBoundaryMapping";
    tests.push_back(test2Daxi_fluxFactor_BoundaryMapping);
    tests.push_back(test3Dtecplot_BoundaryMapping_MSL);
    tests.push_back(test3Dtecplot_BoundaryMapping);
    tests.push_back(test2Daxi_BoundaryMapping);
    tests.push_back(test1D_SurfaceMassEnergyBalance);
    tests.push_back(test1D_BoundaryMapping_constant);
    tests.push_back(test1D_BoundaryMapping_polynomial);
    tests.push_back(test1D_BoundaryMapping_gaussian);
    tests.push_back(test1D_BoundaryMapping_gaussianMixture);
  }

  static TestResult test3Dtecplot_BoundaryMapping_MSL() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: 3D-tecplot";
    std::string testDescription = "Test Boundary Mapping 3D-tecplot";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    list_.append("rhoeUeCH");
    list_.append("h_r");
    word caseName = "test3DMSL";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    label patchID_ = mesh.boundaryMesh().findPatchID("top");

    runTime.setTime(35,0);
    p_.correctBoundaryConditions();
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    scalar expectedScalar_= 208.742612289085741;
    if(!assertEquals((double) p_Bf, (double) expectedScalar_ , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;

    return result;
  }

  static TestResult test3Dtecplot_BoundaryMapping() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: 3D-tecplot";
    std::string testDescription = "Test Boundary Mapping 3D-tecplot";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    list_.append("k");
    word caseName = "test3Dtecplot";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    volTensorField k_
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    runTime.setTime(1,0);
    p_.correctBoundaryConditions();

    scalar p_Bf= p_.boundaryField()[patchID_][0];
    scalar expectedScalar_=10;
    if(!assertEquals((double) p_Bf, (double) expectedScalar_ , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    runTime.setTime(2,0);
    p_.correctBoundaryConditions();

    p_Bf= p_.boundaryField()[patchID_][0];
    expectedScalar_=12782.5486875;

    if(!assertEquals((double) p_Bf,(double) expectedScalar_ , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test2Daxi_BoundaryMapping() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: 2D-axi";
    std::string testDescription = "Test Boundary Mapping 2D-axi";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    list_.append("k");
    word caseName = "test2Daxi";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    volTensorField k_
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    runTime.setTime(0,0);
    p_.correctBoundaryConditions();
    // z = 0.5
    // r = sqrt(1+0.5^2) = 1.1
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    scalar expectedScalar_=15;
    if(!assertEquals(p_Bf, expectedScalar_ , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    runTime.setTime(1,0);
    p_.correctBoundaryConditions();

    p_Bf= p_.boundaryField()[patchID_][0];
    expectedScalar_= 25;
    if(!assertEquals(p_Bf, expectedScalar_ , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test2Daxi_fluxFactor_BoundaryMapping() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: 2Daxi_fluxFactor";
    std::string testDescription = "Test Boundary Mapping 2D-axi fluxFactor";
    TestResult result(suiteName, testName, 1, testDescription);
    word caseName = "test2Daxi_fluxFactor";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField rhoeUeCH_
    (
        IOobject
        (
            "rhoeUeCH",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    runTime.setTime(1,0);

    p_.correctBoundaryConditions();
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    scalar expectedScalar_ = 10123.7791405228;
    if(!assertEquals(p_Bf, expectedScalar_ , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    rhoeUeCH_.correctBoundaryConditions();
    scalar rhoeUeCH_Bf= rhoeUeCH_.boundaryField()[patchID_][0];
    expectedScalar_ = 0.1;
    if(!assertEquals(rhoeUeCH_Bf, expectedScalar_ , &result)) {
      result.expected="rhoeUeCH_.boundaryField()[patchID_][0] = "+name(expectedScalar_);
      result.actual="rhoeUeCH_.boundaryField()[patchID_][0] = "+name(rhoeUeCH_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test1D_SurfaceMassEnergyBalance() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "Bprime";
    std::string testName = "Bprime: constant";
    std::string testDescription = "Test Bprime constant";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    list_.append("k");
    word caseName = "testConstant";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    volScalarField Ta_
    (
        IOobject
        (
            "Ta",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    scalar actual =0;
    scalar expected = 0;
    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    Ta_.correctBoundaryConditions();

    volScalarField& rhoUeCH_ = const_cast<volScalarField&>
                               (
                                   mesh.objectRegistry::lookupObject<volScalarField>("rhoeUeCH")
                               );

    actual = rhoUeCH_.boundaryField()[patchID_][0] ;
    expected = 4.0;

    if(!assertEquals(actual, expected, &result)) {
      result.actual="p_.boundaryField()[patchID_][0] = "+name(actual);
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expected);
      return result;
    }
    runTime.setTime(1,0);
    Ta_.correctBoundaryConditions();
    actual = rhoUeCH_.boundaryField()[patchID_][0] ;
    expected = 5.0;

    if(!assertEquals(actual, expected, &result)) {
      result.actual="p_.boundaryField()[patchID_][0] = "+name(actual);
      result.expected="p_.boundaryField()[patchID_][0] = "+name(expected);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test1D_BoundaryMapping_constant() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: constant";
    std::string testDescription = "Test Boundary Mapping constant";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    list_.append("k");
    word caseName = "testConstant";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    volTensorField k_
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    p_.correctBoundaryConditions();
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 2.0 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 2";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    runTime.setTime(10,0);

    p_.correctBoundaryConditions();
    p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 3.0 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 3";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test1D_BoundaryMapping_polynomial() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: polynomial";
    std::string testDescription = "Test Boundary Mapping polynomial";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    list_.append("k");
    word caseName = "testPolynomial";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    p_.correctBoundaryConditions();
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 1.0 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 1";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }


    runTime.setTime(10,0);

    p_.correctBoundaryConditions();
    p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 21.0 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 21";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test1D_BoundaryMapping_gaussian() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: gaussian";
    std::string testDescription = "Test Boundary Mapping gaussian";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    word caseName = "testGaussian";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    p_.correctBoundaryConditions();
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 4.803947195761616 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 4.803947195761616";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    runTime.setTime(50,0);

    p_.correctBoundaryConditions();
    p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 4.569655926356141 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 4.569655926356141";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }

  static TestResult test1D_BoundaryMapping_gaussianMixture() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "BoundaryMapping";
    std::string testName = "BoundaryMapping: gaussianMixture";
    std::string testDescription = "Test Boundary Mapping gaussian mixture";
    TestResult result(suiteName, testName, 1, testDescription);

    wordList list_;
    list_.append("p");
    word caseName = "testGaussianMixture";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel",caseName);
    word region_ = "porousMat";
    fileName case_ = pato_dir+"/testsuites/MaterialModel/BoundaryConditionsModel/" + caseName;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            region_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    volScalarField p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    label patchID_ = mesh.boundaryMesh().findPatchID("top");
    p_.correctBoundaryConditions();
    scalar p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 10.5500877038 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 10.5500877038";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    runTime.setTime(50,0);

    p_.correctBoundaryConditions();
    p_Bf= p_.boundaryField()[patchID_][0];
    if(!assertEquals(p_Bf, 10.0948400449 , &result)) {
      result.expected="p_.boundaryField()[patchID_][0] = 10.0948400449";
      result.actual="p_.boundaryField()[patchID_][0] = "+name(p_Bf);
      return result;
    }

    Info.level = level_Info_;
    return result;
  }
};

