Info << "\nReading Stratification coefficient"<< endl;
IOdictionary StratProperties
(
 IOobject
 (
  "Stratification",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );

const dimensionedScalar alpha(StratProperties.lookup("alpha"));
