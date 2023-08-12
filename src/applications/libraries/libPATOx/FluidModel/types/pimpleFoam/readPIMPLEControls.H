// We do not have a top-level mesh. Construct the fvSolution for
// the runTime instead.
fvSolution solutionDict(runTime);

const dictionary& pimple = solutionDict.subDict("PIMPLE");

nOuterCorr =
    pimple.lookupOrDefault<int>("nOuterCorrectors", 1);

correctPhi=
    (
        pimple.lookupOrDefault<Switch>("correctPhi", true)
    );
