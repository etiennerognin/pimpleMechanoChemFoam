
// --- Mechanochemical equation
Info<< "Starting mechanochem" << endl;

// Threshold
strength = gradUdotdotA/(epsdot_c*lambda);
// Scission rate (current model)
k = 0.5*pos(strength-1.)*strength*gradUdotdotA;

// Threshold order 0
//k = 10.*pos(strength-1.)/runTime.deltaT();

// Order 1
//k = 1*pos(strength-1.)*gradUdotdotA;

// Micelles model (no threshold)
//k = epsdot_c.value()*gradUdotdotA;


// Degradation equation
fvScalarMatrix cEqn
    (
      fvm::ddt(c)
      + fvm::div(phi,c)
     ==
      - fvm::Sp(k,c)  
    );
cEqn.solve();
c.relax();

// Scission flux based on rate
//kc = k*c;

// Numerical scission flux
kc = - fvc::ddt(c) - fvc::div(phi,c);

// Conversion coefficient
conv = 1. - c;






