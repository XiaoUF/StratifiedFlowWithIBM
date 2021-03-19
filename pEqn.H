volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

MRF.makeRelative(phiHbyA);

surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rho)*mesh.magSf());
phiHbyA += phig;

adjustPhi(phiHbyA, U, p);

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
        p.relax();
        // Correct the momentum source with the pressure gradient flux
        U = HbyA + rAU*fvc::reconstruct((phig - pEqn.flux())/rAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

#include "continuityErrs.H"

