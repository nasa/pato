
#include "PATOx.H"
#include "subtest.h"

#include <vector>


class TestMaterialPropertiesModel : public SubTest
{
 private:
  static void checkRegion(word region, fvMesh& mesh, scalar cp_corr, scalar k_corr, scalar emissivity_corr, scalar absorptivity_corr, TestResult& result) {

    label patchID_ = mesh.boundaryMesh().findPatchID("top");

    simpleEnergyModel& energyModel = meshLookupOrConstructModel<simpleEnergyModel>(mesh,region,"Energy");
    simpleMaterialPropertiesModel& materialPropertiesModel = meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,region,"MaterialProperties");
    materialPropertiesModel.initialize();
    materialPropertiesModel.update();

    scalar cp_Bf = energyModel.refVolField<scalar>("cp").boundaryField()[patchID_][0];
    tensor k_Bf = energyModel.refVolField<tensor>("k").boundaryField()[patchID_][0];
    scalar emissivity_Bf = energyModel.refVolField<scalar>("emissivity").boundaryField()[patchID_][0];
    scalar absorptivity_Bf = energyModel.refVolField<scalar>("absorptivity").boundaryField()[patchID_][0];

    scalar expected_cp = cp_corr;
    scalar expected_k = k_corr;
    scalar expected_emissivity = emissivity_corr;
    scalar expected_absorptivity = absorptivity_corr;

    if(!assertEquals((double) cp_Bf, (double) expected_cp , &result)) {
      result.expected=region+" cp = "+name(expected_cp);
      result.actual="cp = "+name(cp_Bf);
      return;
    }
    if(!assertEquals((double) k_Bf(0,0), (double) expected_k , &result)) {
      result.expected="k(0,0) = "+name(expected_k);
      result.actual="k(0,0) = "+name(k_Bf(0,0));
      return;
    }
    if(!assertEquals((double) k_Bf(1,1), (double) expected_k , &result)) {
      result.expected="k(1,1) = "+name(expected_k);
      result.actual="k(1,1) = "+name(k_Bf(1,1));
      return;
    }
    if(!assertEquals((double) k_Bf(2,2), (double) expected_k , &result)) {
      result.expected="k(2,2) = "+name(expected_k);
      result.actual="k(2,2) = "+name(k_Bf(2,2));
      return;
    }
    if(!assertEquals((double) emissivity_Bf, (double) expected_emissivity , &result)) {
      result.expected="emissivity = "+name(expected_emissivity);
      result.actual="emissivity = "+name(emissivity_Bf);
      return;
    }
    if(!assertEquals((double) absorptivity_Bf, (double) expected_absorptivity , &result)) {
      result.expected="absorptivity = "+name(expected_absorptivity);
      result.actual="absorptivity = "+name(absorptivity_Bf);
      return;
    }

  }
 public:

  TestMaterialPropertiesModel() {
    testSuiteName = "TestMaterialPropertiesModel";
    tests.push_back(test1);
  }

  static TestResult test1() {
    int level_Info_ = Info.level;
    Info.level = 0;
    std::string suiteName = "TestMaterialPropertiesModel";
    std::string testName = "TestMaterialPropertiesModel: Porous Factor";
    std::string testDescription = "Compares material properties obtained using the porous factor module against known values.";
    TestResult result(suiteName, testName, 1, testDescription);

    word caseName = "testPorousFactor";
    const fileName pato_dir = getEnv("PATO_UNIT_TESTING");
    Foam::Time runTime(pato_dir+"/testsuites/MaterialModel/MaterialPropertiesModel",caseName);

    word region_char_ = "char";
    Foam::fvMesh mesh_char
    (
        Foam::IOobject
        (
            region_char_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    word region_virgin_ = "virgin";
    Foam::fvMesh mesh_virgin
    (
        Foam::IOobject
        (
            region_virgin_,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );


    checkRegion(region_char_,mesh_char,855,1.54,0.88,0.8525,result);
    if (!result.failed) {
      checkRegion(region_virgin_,mesh_virgin,810,1.1,0.82125,0.855,result);
    }

    Info.level = level_Info_;

    return result;
  }

};
