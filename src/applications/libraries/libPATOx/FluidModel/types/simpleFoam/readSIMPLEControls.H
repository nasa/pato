// We do not have a top-level mesh. Construct the fvSolution for
// the runTime instead.
fvSolution solutionDict(runTime);

//const dictionary& simple = solutionDict.subDict("SIMPLE");

//correctPhi=        (
//    simple.lookupOrDefault<Switch>("correctPhi", true)
//                   );
