CoNum = -GREAT;

forAll(fluidRegions, regionI)
{
  CoNum =
      max
      (
          courantNo
          (
              fluidRegions[regionI],
              runTime,
              phiFluid[regionI]
          ),
          CoNum
      );
}
