/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  temperatureDependentViscosity;

nu table
(
    (869 5.87866E-07)
    (1501.67 3.52427E-07)
    (2134.33 3.01095E-07)
    (2767 2.89045E-07)
);

rho 2409.86;

kappa table
(
    (869.00 65.86)
    (2767.00 105.24)
);

kappas table
(  
    (298.10 162.05)
    (375.99 165.71)
    (830.00 135.26)
);


cp table
(   
    (869.00 1188.00)
    (2767.00 928.94)
);

cps table
(   
    (297.75 977.95)
    (830.00 1195.24)
);
    
Tsolidus 806.0;
Tliquidus 925.0;
LatentHeat 494000;
beta 0.000025;

elec_resistivity 2.40e-7;

// ************************************************************************* //
