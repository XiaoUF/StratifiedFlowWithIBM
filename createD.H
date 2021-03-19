
// Read Darcy's coefficient
Info << "\nReading Darcy's coefficient"<< endl;
IOdictionary DarcyProperties
(
 IOobject
 (
  "DarcyProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );

// K is the permeability coefficient for the momentum equation
const dimensionedScalar K(DarcyProperties.lookup("K"));
// xc is the center of the hill
const dimensionedScalar xc(DarcyProperties.lookup("xc"));
// Hi is the height of the hill
const dimensionedScalar Hi(DarcyProperties.lookup("Hi"));
// Li is the length for the Gaussian shape
const dimensionedScalar Li(DarcyProperties.lookup("Li"));
// L0i is the length for the cutoff
const dimensionedScalar L0i(DarcyProperties.lookup("L0i"));
// ks is the coefficient for the Salinity equation
const dimensionedScalar KS(DarcyProperties.lookup("KS"));
// U0 is the moving/towing speed
const dimensionedScalar U0(DarcyProperties.lookup("U0"));
// Sponge layer parameter
const dimensionedScalar HsL(DarcyProperties.lookup("HsL"));
const dimensionedScalar HsH(DarcyProperties.lookup("HsH"));
const dimensionedScalar DSL(DarcyProperties.lookup("DSL"));
//const dimensionedVector f_cur(DarcyProperties.lookup("fc"));
// This is used to create the coefficients required for Porous media

volScalarField D
(
 IOobject
 (
  "D_pf",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("0",K.dimensions(),0.0)
 );

volScalarField DS
(
 IOobject
 (
  "D_Sal",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("0",KS.dimensions(),0.0)
 );

volVectorField Uf
(
 IOobject
 (
  "Uf",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedVector("0",U.dimensions(),vector(0.0,0.0,0.0))
 );

volScalarField Sf
(
 IOobject
 (
  "Sf",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("0",S.dimensions(),0.0)
 );

//create the permeability 
Info << "Setting the permeability...\n" << endl;
forAll(mesh.C(),i)
{
  vector coords = mesh.C()[i];
  scalar Hz;
  scalar l2;
  scalar xc0;
  xc0 = xc.value()+U0.value()*runTime.value();
  l2 = -0.5*(coords[0]-xc0)*(coords[0]-xc0);
  l2 = l2/Li.value()/Li.value();
  Hz = Hi.value()*Foam::exp(l2);
  Uf[i] = vector(0.0,0.0,0.0);
  D[i]  = 0.0;
  DS[i] = 0.0;
  if(coords[2]<Hz){
    Uf[i] = vector(U0.value(),0,0);
    D[i]  = K.value();
    DS[i] = KS.value();
  }
  scalar h1 = HsL.value();
  scalar h2 = HsH.value();
  if(coords[2]>h1){
    D[i] = DSL.value()*(coords[2]-h1)/(h2-h1);
  }
//  Info << "rsq = " << rsq << "; dist = " << dist << endl;
//  if(dist <= rsq){
//    Info << "dist = " << dist << endl;
//    D[i] = K.value();
//  }
}
