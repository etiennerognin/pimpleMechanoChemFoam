
// By looping inside the outer PIMPLE loop, the elastic stress is 
// coupled with the velocity field.

// The dev function is used to ensure that no artificial expansion is
// introduced in the conformation equation
gradU = dev(fvc::grad(U, "D")); 

// For visualisation only
D = symm(gradU);

gradUdotdotA = (gradU && A);
// Flooring at zero
dimensionedScalar zeroT("zeroT",dimensionSet(0,0,-1,0,0,0,0),0.);
gradUdotdotA = max(gradUdotdotA, zeroT);
            
// Equation for A
//
dimensionedSymmTensor Id("Id",dimensionSet(0,0,0,0,0,0,0),symmTensor::I);

fvSymmTensorMatrix AEqn
(
    fvm::ddt(A)  
  + fvm::div(phi, A)
 ==                
  twoSymm(A & gradU)              
  // Finite extensibility model
  - fvm::Sp(2.*(1.75-0.75*tr(A))*gradUdotdotA, A)                
  // Zimm relaxation              
  - fvm::Sp(1./tau_Z, A)  
  // Equilibirum value                     
  + 1./(3*sqr(xi)*tau_Z)*Id                               
);

// Evolve A
AEqn.solve();

// Correct A (better ways to keep positiveness exist)
dimensionedSymmTensor Amin("Amin",dimensionSet(0,0,0,0,0,0,0),symmTensor(0,-1,-1,0,-1,0));
A = max(A, Amin);

// update gradU:A
gradUdotdotA = max(gradU && A, zeroT);
lambda = max(tr(A), 1./sqr(xi)); 
            
