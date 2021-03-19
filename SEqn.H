// Scalar transport equation, simple turbulent advection-diffusion equation
{
fvScalarMatrix SEqn
  (
   fvm::ddt(S)
   +fvm::div(phi,S)
   -fvm::laplacian(turbulence->nuEff(),S)
   + fvm::Sp(DS,S)
   );

 SEqn.relax();
 SEqn.solve();

 S.correctBoundaryConditions();
}
