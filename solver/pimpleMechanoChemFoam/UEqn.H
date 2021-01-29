// Solve the Momentum equation

// Stress-chemistry mixture model
volScalarField eta_p_coef = ((1.6188844-1.)*c+1.)/1.6188844;
volScalarField eta_E_coef = (3.*c+1.)/4.;

// Polymeric stress tensor
tau_p = (3*sqr(xi)*eta_p_coef*eta_p/tau_Z + eta_E_coef*eta_E*gradUdotdotA)*A;

// Polymer force
f_p = fvc::div(tau_p);

// Effective viscosity (for record only)
dimensionedScalar SMALLD2("SMALLD2",dimensionSet(0,0,-2,0,0,0,0),1e-9);
eta_eff = 0.5*(3*sqr(xi)*eta_p/tau_Z + eta_E*gradUdotdotA)*(gradU && A)/((D && gradU) + SMALLD2);

// -- BSD technique --
// Equivalent surface viscosity
surfaceScalarField eta_starf = eta_E*((mesh.Sf() & fvc::interpolate(lambda*A)) & mesh.Sf())/ magSqr(mesh.Sf());

// Avoid negative viscosity (occurence means A is not semi-positive anymore)
eta_starf = mag(eta_starf);

// For paraFoam plotting
eta_star1 = fvc::average(eta_starf);

// --- Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + MRF.DDt(U)
  // Solvent viscosity
  - fvm::laplacian(eta_s/rho, U)

  // Polymer explicit contribution
  - f_p/rho
  
  // Both Sides Diffusion  
  - fvm::laplacian(eta_starf/rho, U, "BSD")
  + fvc::surfaceIntegrate(eta_starf/rho*mesh.Sf() & fvc::interpolate(fvc::grad(U, "BSD"))) 
);                                       

UEqn.relax();

if (pimple.momentumPredictor() or firstIter)
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
