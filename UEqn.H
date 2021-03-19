// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  + fvm::Sp(D,U)
  - D*Uf
 ==
    fvOptions(U)
    
);

UEqn.relax();

rho = 1. + alpha*S;

fvOptions.constrain(UEqn);

{
    solve(UEqn == 
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p)
                )*mesh.magSf()
            )
         );

    fvOptions.correct(U);
}

