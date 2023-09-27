const dictionary& stressControl =
    mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
Switch predictor(stressControl.lookup("predictor"));
int infoFrequency(readInt(stressControl.lookup("infoFrequency")));
Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
