#include "createRDeltaT.H"

const dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);

rho.oldTime();


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p +=
        dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureReference.refValue()
          - getRefCellValue(p, pressureReference.refCell())
        );

    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

volScalarField& alpha1(mixture.alpha1());

mesh.schemes().setFluxRequired(alpha1.name());


IOdictionary transportProperties_metal
(
    IOobject
    (
        "physicalProperties.metal",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary transportProperties_gas
(
    IOobject
    (
        "physicalProperties.gas",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar Tsolidus1
(
    "Tsolidus1",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties_metal.lookup("Tsolidus")
);

const dimensionedScalar Tsolidus2
(
    "Tsolidus2",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties_gas.lookup("Tsolidus")
);

const dimensionedScalar Tliquidus1
(
    "Tliquidus1",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties_metal.lookup("Tliquidus")
);

const dimensionedScalar Tliquidus2
(
    "Tliquidus2",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties_gas.lookup("Tliquidus")
);

const dimensionedScalar LatentHeat1
(
    "LatentHeat1",
    dimensionSet(0, 2, -2, 0, 0),
    transportProperties_metal.lookup("LatentHeat")
);

const dimensionedScalar LatentHeat2
(
    "LatentHeat2",
    dimensionSet(0, 2, -2, 0, 0),
    transportProperties_gas.lookup("LatentHeat")
);

const dimensionedScalar beta1
(
    "beta1",
    dimensionSet(0, 0, 0, -1, 0),
    transportProperties_metal.lookup("beta")
);

const dimensionedScalar beta2
(
    "beta2",
    dimensionSet(0, 0, 0, -1, 0),
    transportProperties_gas.lookup("beta")
);

const dimensionedScalar Marangoni_Constant
(
    "Marangoni_Constant",
    dimensionSet(1, 0, -2, -1, 0),
    phaseProperties.lookup("dsigmadT")
);

const dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    phaseProperties.lookup("p0")
); // atmospheric pressure

const dimensionedScalar Tvap
(
    "Tvap",
    dimensionSet(0, 0, 0, 1, 0),
    phaseProperties.lookup("Tvap")
);

const dimensionedScalar Mm
(
    "Mm",
    dimensionSet(1, 0, 0, 0, -1),
    phaseProperties.lookup("Mm")
);

const dimensionedScalar LatentHeatVap
(
    "LatentHeatVap",
    dimensionSet(0, 2, -2, 0, 0),
    phaseProperties.lookup("LatentHeatVap")
);

const dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.314
);

/*const dimensionedScalar kappa1l
(
    "kappa1l",
    dimensionSet(1, 1, -3, -1, 0),
    transportProperties_metal.lookup("kappa")
);

const dimensionedScalar kappa1s
(
    "kappa1s",
    dimensionSet(1, 1, -3, -1, 0),
    transportProperties_metal.lookup("kappasolid")
);

const dimensionedScalar Cp1l
(
    "Cp1l",
    dimensionSet(0, 2, -2, -1, 0),
    transportProperties_metal.lookup("Cp")
); */

//const Polynomial<8> polykappa1l(transportProperties_metal.lookup("poly_kappa"));
//const Polynomial<8> polykappa1s(transportProperties_metal.lookup("poly_kappasolid"));
//const Polynomial<8> polycp1l(transportProperties_metal.lookup("poly_Cp"));
//const Polynomial<8> polycp1s(transportProperties_metal.lookup("poly_Cpsolid"));

autoPtr<Function1<scalar>> kappaPtr = Function1<scalar>::New("kappa", transportProperties_metal);
autoPtr<Function1<scalar>> kappasPtr = Function1<scalar>::New("kappas", transportProperties_metal);
autoPtr<Function1<scalar>> cpPtr = Function1<scalar>::New("cp", transportProperties_metal);
autoPtr<Function1<scalar>> cpsPtr = Function1<scalar>::New("cps", transportProperties_metal);
//Function1<scalar> kappa_Ptr("kappa_list", transportProperties_metal);


/*
const dimensionedScalar Cp1s
(
    "Cp1s",
    dimensionSet(0, 2, -2, -1, 0),
    transportProperties_metal.lookup("Cpsolid")
);*/

const dimensionedScalar kappa2l
(
    "kappa2l",
    dimensionSet(1, 1, -3, -1, 0),
    transportProperties_gas.lookup("kappa")
);

const dimensionedScalar kappa2s
(
    "kappa2s",
    dimensionSet(1, 1, -3, -1, 0),
    transportProperties_gas.lookup("kappasolid")
);

const dimensionedScalar Cp2l
(
    "Cp2l",
    dimensionSet(0, 2, -2, -1, 0),
    transportProperties_gas.lookup("Cp")
);

const dimensionedScalar Cp2s
(
    "Cp2s",
    dimensionSet(0, 2, -2, -1, 0),
    transportProperties_gas.lookup("Cpsolid")
);

//const Polynomial<8> polykappa_m(transportProperties_metal.lookup("poly_kappa"));
//const Polynomial<8> polykappa_g(transportProperties_gas.lookup("poly_kappa"));

//const Polynomial<8> polycp_m(transportProperties_metal.lookup("poly_cp"));
//const Polynomial<8> polycp_g(transportProperties_gas.lookup("poly_cp"));


volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
);

//#include "updateKappaCp.H"

volScalarField TSolidus
(
    IOobject
    (
        "TSolidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2
);

volScalarField TLiquidus
(
    IOobject
    (
        "TLiquidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2
);

volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DC", dimensionSet(1,-3,-1,0,0), 1.0e14)
);

volScalarField metalaverage
(
    IOobject
    (
        "metalaverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::average(alpha1)
);

volScalarField epsilon1
(
    IOobject
    (
        "epsilon1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(min((T - TSolidus)/(TLiquidus - TSolidus), 1.0), 0.0)
);

#include "updateKappaCp.H"

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered", dimless, vector::zero)
);

volScalarField epsilon1mask
(
    IOobject
    (
        "epsilon1mask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask",dimless,0.0)
);

volScalarField epsilon1mask2
(
    IOobject
    (
        "epsilon1mask2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask2", dimless, 0.0)
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0), vector::zero)
);

volScalarField LatentHeat
(
    IOobject
    (
        "LatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*LatentHeat1 + (1.0 - alpha1)*LatentHeat2
);

volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

volScalarField Tcorr
(
    IOobject
    (
        "Tcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tcorr", dimensionSet(0, 0, 0, 1, 0), 0.0)
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*beta1 + (1.0-alpha1)*beta2
);

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhok", dimless, 0.0)
);

volScalarField Num_divU
(
    IOobject
    (
        "Num_divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Num_divU", dimensionSet(0,0,-1,0,0), 0.0)
);

volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

// Recoil Pressure
volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

// Evaporation Cooling
volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);

// volScalarField ddte1
// (
//     IOobject
//     (
//         "ddte1",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("ddte1",dimensionSet(0,0,-1,0,0),0.0),
//     zeroGradientFvPatchScalarField::typeName
// );

volScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("damper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDamper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volVectorField gradepsilon1(fvc::grad(alpha_filtered));
volScalarField e1temp(fvc::average(epsilon1));
volScalarField e2temp(fvc::average(epsilon1mask));

dimensionedScalar elec_resistivity_metal
(
    "elec_resistivity",
    dimensionSet(1, 3, -3, 0, 0),
    transportProperties_metal.lookup("elec_resistivity")
);

// Electrical resistivity of condensate only
volScalarField electrical_resistivity
(
    IOobject
    (
        "electrical_resistivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    elec_resistivity_metal
);

#include "createAlphaFluxes.H"

// Construct incompressible turbulence model
incompressibleInterPhaseTransportModel turbulence
(
    U,
    phi,
    alphaPhi1,
    mixture
);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// Create laser heat source
laserHeatSource laser(mesh);

volScalarField rhoCp(rho*cp);
rhoCp.oldTime();

surfaceScalarField rhoPhi("rhoPhi", fvc::interpolate(rho)*phi);
surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

#include "updateProps.H"

const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless, 1.0e-12
);

volVectorField gradAlpha(fvc::grad(alpha1));

volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));

// This variable is >=0 if the cell has been liquid
volScalarField meltHistory
(
    IOobject
    (
        "meltHistory",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("meltHistory", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

// This variable is 1 if (alpha.metal >= 0.5 and epsilon1 >=0.5), 0 otherwise
volScalarField condition
(
    IOobject
    (
        "condition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("condition", dimensionSet(0, 0, 0, -0, 0), 0.0)
);


// added by MT
volScalarField depth_y
(
    IOobject
    (
        "depth_y",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),      
    mesh,       
    dimensionedScalar("depth_y", dimensionSet(0, 0, 0, 0, 0), 0.0)
);                  
                    
volScalarField melt_depth_y
(                   
    IOobject        
    (             
        "melt_depth_y",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),  
    mesh,
    dimensionedScalar("melt_depth_y", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField magU
(             
    IOobject        
    (        
        "magU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),  
    mesh,
	dimensionedScalar("magU", dimensionSet(0, 1, -1, 0, 0), 0.0)
);

// ### added evaporation part by MT
// 1. Create X for each metal element Al Cu Mg Zn
int metal_comp_num = 4;
List<word> element_name(metal_comp_num);
element_name[0] = "Al";
element_name[1] = "Cu";
element_name[2] = "Mg";
element_name[3] = "Zn";


PtrList<volScalarField> Xs(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    Xs.set
    (
        i, 
        new volScalarField
        (   
            IOobject
            (
                "X."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
                mesh,
                dimensionedScalar("X", dimensionSet(0, 0, 0, 0, 0), 0.0)
        )   
    );

}



Info << "Xs created. ( ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << Xs[i].name() << " ";
}
Info << ")" << endl;

PtrList<volScalarField> Ws(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    Ws.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "W."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
                mesh,
                dimensionedScalar("W", dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    );

}
        

        
Info << "Ws created. ( ";
for (int i = 0; i < metal_comp_num; i++)
{   
    Info << Ws[i].name() << " ";
}       
Info << ")" << endl;


// 2. Create variables
double* total_evap;
total_evap = new double[metal_comp_num];

for (int i = 0; i < metal_comp_num; i++)
{
    total_evap[i] = 0.0;
}

double* total_evap2;
total_evap2 = new double[metal_comp_num];
        
for (int i = 0; i < metal_comp_num; i++)
{       
    total_evap2[i] = 0.0;
}

double* total_evap3;
total_evap3 = new double[metal_comp_num];

for (int i = 0; i < metal_comp_num; i++)
{
    total_evap3[i] = 0.0;
}

scalar total_evap_air = 0.0;


// activity coefficient

PtrList<volScalarField> gammas(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    gammas.set
    (   
        i,
        new volScalarField
        (
            IOobject
            (
                "gamma."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("gamma", dimensionSet(0, 0, 0, 0, 0), 1.00)
        )
    );
}


Info << "gammas is created. (";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << gammas[i].name() << " ";
}
Info << ")" << endl;

// Set evaporation related terms 

volScalarField molarmass_total
(
    IOobject
    (
        "molarmass_total",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("molarmass_total", dimensionSet(1, 0, 0, 0, -1), 0.0)
);

volScalarField Ps_total
(
    IOobject
    (
        "Ps.total",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ps.total", dimensionSet(1, -1, -2, 0, 0), 0.0)
);

PtrList<volScalarField> Pss(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    Pss.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Ps."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ps", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );
}



Info << "Pss created. ( ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << Pss[i].name() << " ";
}
Info << ")" << endl;



volScalarField M3_total
(
    IOobject
    (
        "M3.total",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("M3.total", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Phi", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField Alpha
(
    IOobject
    (
        "Alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Alpha", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

// alphas and volume fraction of metals for changes in composition
PtrList<volScalarField> alphas(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    alphas.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alphas."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("alpha", dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    );
}

for (int i = 0; i < metal_comp_num; i++)
{
    forAll(alphas[i].boundaryField(), patchi)
    {
        alphas[i].boundaryFieldRef().set
        (
            patchi,
            fvPatchField<scalar>::New
            (
                "zeroGradient",
                mesh.boundary()[patchi],
                alphas[i]
            )
        );
    }
}


Info << "alphas created. ( ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << alphas[i].name() << " ";
}
Info << ")" << endl;


PtrList<volScalarField> Cs(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    Cs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "C."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("C", dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    );
}

for (int i = 0; i < metal_comp_num; i++)
{
    forAll(Cs[i].boundaryField(), patchi)
    {
        Cs[i].boundaryFieldRef().set
        (
            patchi,
            fvPatchField<scalar>::New
            (
                "zeroGradient",
                mesh.boundary()[patchi],
                Cs[i]
            )
        );
    }
}



Info << "Cs created. ( ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << Cs[i].name() << " ";
}
Info << ")" << endl;

volScalarField total_mol_vol
(               
    IOobject
    (
        "total_mol_vol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,       
    dimensionedScalar("total_mol_vol", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField total_weight
(
    IOobject
    (
        "total_weight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("total_weight", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

int n_comp = metal_comp_num;

volScalarField zero_vol
(
    IOobject
    (
        "zero_vol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero_vol", dimensionSet(0, 0, -1, 0, 0), 0.0)
);


volScalarField temp_vol
(
    IOobject
    (
        "temp_vol", 
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp_vol", dimensionSet(0, 0, -1, 0, 0), -1e-5)
);

PtrList<volScalarField> massfluxs(metal_comp_num);

for (int i = 0; i < metal_comp_num; i++)
{
    massfluxs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "massflux."+element_name[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("massfluxC", dimensionSet(0, 1, -1, 0, 0), 0.0)
        )
    );
}


Info << "massfluxs created. ( ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << massfluxs[i].name() << " ";
}
Info << ")" << endl;


volScalarField total_massflux
(
	IOobject
	(
		"total_massflux",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("massfluxC", dimensionSet(0, 1, -1, 0, 0), 0.0)
);



volScalarField melted
(
    IOobject
    (
        "melted",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("melted", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField melted_metal
(   
    IOobject
    (   
        "melted_metal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("melted_metal", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

// 3. L parameter for binary and ternalry (Al Ti V)
// must fix this part
// binary
double FL1[6][4][4];
for (int i1 = 0; i1 < 6; i1++)
    for (int i2 = 0; i2 < 4; i2++)
        for (int i3 = 0; i3 < 4; i3++)
        {
            FL1[i1][i2][i3] = 0.0;
        }

// Al Cu Mg Zn
// Al Cu
FL1[0][0][1] = -67094;
FL1[1][0][1] = 8.555;
FL1[2][0][1] = 32148;
FL1[3][0][1] = -7.118;
FL1[4][0][1] = 5915;
FL1[5][0][1] = -5.889;

// Al Mg
FL1[0][0][2] = -12000;
FL1[1][0][2] = 8.566;
FL1[2][0][2] = 1894;
FL1[3][0][2] = -3;
FL1[4][0][2] = 2000;

// Al Zn
FL1[0][0][3] = 10465.55;
FL1[1][0][3] = -3.39259;

// Cu Mg
FL1[0][1][2] = -36984;
FL1[1][1][2] = 4.7561;
FL1[2][1][2] = -8191.29;

// Cu Zn
FL1[0][1][3] = -40688.6;
FL1[1][1][3] = 12.8655;
FL1[2][1][3] = 2252.3;
FL1[3][1][3] = -4.8069;
FL1[4][1][3] = 5405.9;
FL1[5][1][3] = -2.9197;

// Mg Zn
FL1[0][2][3] = -77729.24;
FL1[1][2][3] = 680.52266;
FL1[2][2][3] = 3674.72;
FL1[3][2][3] = 0.57139;
FL1[4][2][3] = -1588.15;


// ternary
double TL1[8][4][4][4];
for (int i1 = 0; i1 < 8; i1++)
    for (int i2 = 0; i2 < 4; i2++)
        for (int i3 = 0; i3 < 4; i3++)
            for (int i4 = 0; i4 < 4; i4++)
            {
                TL1[i1][i2][i3][i4] = 0.0;
            }

// Al Cu Mg Zn
    // TL1[8][n_comp][n_comp][n_comp]  // adding temperature dependence of ternary parameter
// Al Cu Mg - X

// Al Cu Zn
TL1[0][0][1][3] = 0;
TL1[1][0][1][3] = 0;
TL1[2][0][1][3] = 5475;
TL1[3][0][1][3] = 13000;
TL1[4][0][1][3] = 43800;
TL1[5][0][1][3] = 0;
TL1[6][0][1][3] = 0;
TL1[7][0][1][3] = 0;

// Al Mg Zn
TL1[0][0][2][3] = 0;
TL1[1][0][2][3] = 0;
TL1[2][0][2][3] = -11475;
TL1[3][0][2][3] = -11475;
TL1[4][0][2][3] = -11475;
TL1[5][0][2][3] = 11;
TL1[6][0][2][3] = 11;
TL1[7][0][2][3] = 11;

// Cu Mg Zn - X



// 4. properties
// must fix this part
double Latentheatgas_list[4] = {9462000.0, 4720995.0, 5241719.8, 1754000};
double Tvap_list[4] = {2790, 2835, 1363, 1180};
double molarmass_list[4] = {0.02698, 0.063546, 0.024305, 0.06538};
double rho_list[4] = {2698.9, 8000, 1590, 6580};
double wf_list[4] = {0.9015, 0.0152, 0.0244, 0.0589};
scalar air_molarmass = 0.039948;

// no fix
double xf_list[4] = {0.0, 0.0, 0.0, 0.0};
double vf_list[4] = {0.0, 0.0, 0.0, 0.0};


// 5. Set compositions and print (wf -> mf(X) -> vf(C))
scalar total_mole = 0.0;
scalar total_vol = 0.0;

for (int i = 0; i < metal_comp_num; i++)
{
	// total mole
	total_mole += wf_list[i] / molarmass_list[i];
    // total vol
    total_vol += wf_list[i] / rho_list[i];
}

for (int i = 0; i < metal_comp_num; i++)
{
    Xs[i] = wf_list[i] / molarmass_list[i] / total_mole;
    Cs[i] = wf_list[i] / rho_list[i] / total_vol;
    xf_list[i] = Xs[i][0];
    vf_list[i] = Cs[i][0];
}

for (int i = 0; i < metal_comp_num; i++)
{
    Info << "  X[" << i << "] = " << Xs[i][0];
}
Info << endl;

for (int i = 0; i < metal_comp_num; i++)
{
    Info << "  C[" << i << "] = " << Cs[i][0];
}
Info << endl;

Info << "mole fraction reference: ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << "  xf[" << i << "] = " << xf_list[i];
}
Info << endl;

Info << "volume fraction reference: ";
for (int i = 0; i < metal_comp_num; i++)
{
    Info << "  vf[" << i << "] = " << vf_list[i];
}
Info << endl;

Info << " " << endl;

// set_alphas
for (int i = 0; i < metal_comp_num; i++)
{
    alphas[i] = alpha1 * vf_list[i];
}

