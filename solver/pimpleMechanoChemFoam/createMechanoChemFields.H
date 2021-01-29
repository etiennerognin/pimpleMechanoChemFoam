// Solved fields ---------------------------------------------------------------

// Normalized Conformation tensor
volSymmTensorField A
(
    IOobject
    (
        "A",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Concentration of intact chains
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Physicochemical properties --------------------------------------------------

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Solvent viscosity
dimensionedScalar eta_s
(
    "eta_s",
    dimensionSet(1,-1,-1,0,0,0,0),
    transportProperties.lookup("eta_s")
);

// Fluid density
dimensionedScalar rho
(
    "rho",
    dimensionSet(1,-3,0,0,0,0,0),
    transportProperties.lookup("rho")
);


// Relaxation time (Zimm time)
dimensionedScalar tau_Z
(
    "tau_Z",
    dimensionSet(0,0,1,0,0,0,0),
    transportProperties.lookup("tau_Z")
);

// Extensibility of polymer
dimensionedScalar xi
(
    "xi",
    dimensionSet(0,0,0,0,0,0,0),
    transportProperties.lookup("xi")
);

// Added low-shear viscosity
dimensionedScalar eta_p
(
    "eta_p",
    dimensionSet(1,-1,-1,0,0,0,0),
    transportProperties.lookup("eta_p")
);

// Added extensional viscosity
dimensionedScalar eta_E
(
    "eta_E",
    dimensionSet(1,-1,-1,0,0,0,0),
    transportProperties.lookup("eta_E")
);

// Critical strain rate epsdot_c
dimensionedScalar epsdot_c
(
    "epsdot_c",
    dimensionSet(0,0,-1,0,0,0,0),
    transportProperties.lookup("epsdot_c")
);

// Calculated fields -----------------------------------------------------------
// Velocity gradient
volTensorField gradU
(
    IOobject
    (
        "gradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("gradU",dimensionSet(0,0,-1,0,0,0,0),tensor::zero)
);


// Strain rate tensor
volSymmTensorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("D",dimensionSet(0,0,-1,0,0,0,0),symmTensor::zero)
);


// Polymer stress
volSymmTensorField tau_p
(
    IOobject
    (
        "tau_p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("tau_p",dimensionSet(1,-1,-2,0,0,0,0),symmTensor::zero)
);

volVectorField f_p
(
    IOobject
    (
        "f_p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("f_p",dimensionSet(1,-2,-2,0,0,0,0),vector::zero)
);

// Polymer extension
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda",dimensionSet(0,0,0,0,0,0,0),0.)
);

// Polymer friction
volScalarField gradUdotdotA
(
    IOobject
    (
        "gradUdotdotA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gradUdotdotA",dimensionSet(0,0,-1,0,0,0,0),0.)
);

// Polymer effective viscosity
volScalarField eta_eff
(
    IOobject
    (
        "eta_eff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("eta_eff",dimensionSet(1,-1,-1,0,0,0,0),0.)
);

// BSD trace of viscosity
volScalarField eta_star1
(
    IOobject
    (
        "eta_star1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("eta_star1",dimensionSet(1,-1,-1,0,0,0,0),0.)
);

// First order scission rate
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("k",dimensionSet(0,0,-1,0,0,0,0),0.)
);

// Conversion factor
volScalarField conv
(
    IOobject
    (
        "conv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("conv",dimensionSet(0,0,0,0,0,0,0),0.)
);

// Conversion rate
volScalarField kc
(
    IOobject
    (
        "kc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kc",dimensionSet(0,0,-1,0,0,0,0),0.)
);

// Elongation rate along polymer principal direction divided by critical rate
volScalarField strength
(
    IOobject
    (
        "strength",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("strength",dimensionSet(0,0,0,0,0,0,0),0.)
);

