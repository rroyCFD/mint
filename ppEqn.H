{
    // The inverse of the momentum equation "A" operator (the matrix diagonal) and
    // its projection from cell centers to cell faces.  Because it is a diagonal matrix,
    // The inverse is just the diagonal matrix of the reciprocals, hence the "r" in the
    // name "rAU".
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    // The momentum equation "H" operator is the off-diagonal times the last known U.
    // In the equations it is always multiplied with inv(A).  This is HbyA.
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    // This is the Boussinesq buoyancy term multipled by the inverse of the A operator.
    Boussinesq.updateBuoyancyTerm();
    surfaceScalarField phig(rAUf * Boussinesq.buoyancyTerm() * mesh.magSf());

    // Project HbyA to cell faces and apply a correction for time stepping scheme.
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U, phi)
    );

    // This part balances global mass flux.  It does it in a temporary field, and then
    // applies the correction indirectly by setting the pressure gradient to be used in
    // the fixedFluxPressure boundary condition on p_rgh, or directly if the zeroGradient
    // boundary condition on p_rgh is used.
    surfaceScalarField phiFixedFlux = phiHbyA;

    adjustPhi(phiFixedFlux, U, p_rgh);

    forAll(p_rgh.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p_rgh.boundaryField()[patchi]))
        {
            phiHbyA.boundaryFieldRef()[patchi] = phiFixedFlux.boundaryField()[patchi];
        }
    }

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - phiFixedFlux.boundaryField()
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );
    p_rgh.correctBoundaryConditions();

    // Apply the same gradient to pp following the coresponding p_rgh BCs
    forAll(p_rgh.boundaryField(), patchi)
    {
        if(isType<fixedFluxPressureFvPatchScalarField> (p_rgh.boundaryField()[patchi]))
        {
            fixedFluxPressureFvPatchScalarField& p_rghGradPatch = refCast<fixedFluxPressureFvPatchScalarField>(p_rgh.boundaryFieldRef()[patchi]);
            fixedGradientFvPatchScalarField& ppGradPatch        = refCast<fixedGradientFvPatchScalarField>(pp.boundaryFieldRef()[patchi]);

            ppGradPatch.gradient() = p_rghGradPatch.gradient();
            Info << "minGrad: " << gMin(ppGradPatch.gradient()) << tab
                 << "maxGrad: " << gMax(ppGradPatch.gradient()) << endl;
        }
    }

    phipp = phiHbyA - (rAUf * fvc::snGrad(p_rgh) * mesh.magSf());


    // Non-orthogonal corrector loop.
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix ppEqn
        (
            fvm::laplacian(rAUf, pp) == fvc::div(phipp)
        );

        ppEqn.setReference(pRefCell, 0.0);

        ppEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phipp -= ppEqn.flux();

            // Explicitly relax pressure for momentum corrector
            pp.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            Upp = U -rAU*fvc::grad(pp);
            Upp.correctBoundaryConditions();
            fvOptions.correct(Upp);
        }
    }

    // #include "adjustPressureLevel.H"

    volVectorField Udiff("Udiff", U - Upp);
    Info << Udiff.name() << "\t min: " << gMin(Udiff) << "\tmax: " << gMax(Udiff) << endl;
}
