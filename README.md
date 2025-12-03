# EXP 3 B: IIR-CHEBYSHEV-FITER-DESIGN

## AIM: 

 To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 
```c
clc;
clear;
close;

// ------------------------------------
// Chebyshev Type-I IIR Low-Pass Filter
// ------------------------------------

// ---- Filter Specifications (EDIT HERE) ----
wp = 0.5*%pi;      // passband frequency (rad)
ws = 0.8*%pi;      // stopband frequency (rad)
alphap = 4;        // passband ripple (dB)
alphas = 20;       // stopband attenuation (dB)
T = 1;             // sampling time

// ---- Prewarp ----
omegap = (2/T) * tan(wp/2);
omegas = (2/T) * tan(ws/2);

// ---- Filter Order ----
N = acosh( sqrt((10^(0.1*alphas)-1)/(10^(0.1*alphap)-1)) ) / acosh(omegas/omegap);
N = ceil(N);

// ---- Epsilon & Cutoff ----
Eps = sqrt(10^(0.1*alphap) - 1);
omegac = omegap / (Eps)^(1/N);

// ---- Analog Chebyshev LPF ----
[pols, gn] = zpch1(N, Eps, omegac);
hs = poly(gn, 's', 'coeff') / real(poly(pols, 's'));

// ---- Bilinear Transform H(z) ----
z = poly(0,'z');
Hz = horner(hs, (2/T) * ((z-1)/(z+1)));

// ---- Frequency Response ----
Nf = 512;
HW = frmag(Hz, Nf);
w = 0 : %pi/(Nf-1) : %pi;

// ---- Plot ----
figure;
plot(w/%pi, abs(HW));
xlabel("Normalized Frequency (×π rad/sample)");
ylabel("Magnitude");
title("Chebyshev Type-I IIR Low-Pass Filter Frequency Response");
xgrid();
```


## PROGRAM (HPF): 
```c
clc;
clear;
close;

// ------------------------------------
// Chebyshev Type-I IIR High-Pass Filter
// ------------------------------------

// ----- Filter Specifications (EDIT THESE) -----
wp = 0.8*%pi;      // passband edge (rad)
ws = 0.5*%pi;      // stopband edge (rad)
alphap = 4;        // passband ripple (dB)
alphas = 20;       // stopband attenuation (dB)
T = 1;             // sampling time

// ----- Prewarp -----
Omegap = (2/T) * tan(wp/2);
Omegas = (2/T) * tan(ws/2);

// ----- For HPF: Convert to equivalent LPF specs -----
Omegar = Omegap;                      // reference for HPF
Omegap_lpf = Omegar / Omegap;
Omegas_lpf = Omegar / Omegas;

// ----- Filter Order -----
N = acosh( sqrt((10^(0.1*alphas)-1)/(10^(0.1*alphap)-1)) ) / acosh(Omegas_lpf / Omegap_lpf);
N = ceil(N);

// ----- Epsilon & Cutoff -----
Eps = sqrt(10^(0.1*alphap) - 1);
Omegac_lpf = Omegap_lpf / (Eps)^(1/N);

// ----- Analog Low-Pass Prototype -----
[pols, gn] = zpch1(N, Eps, Omegac_lpf);
Hs_lpf = poly(gn, "s", "coeff") / real(poly(pols, "s"));

// ----- Low-Pass → High-Pass Transformation -----
s = poly(0, "s");
H_hp_s = horner(Hs_lpf, Omegac_lpf ./ s);   // s -> Ωc/s substitution

// ----- Bilinear Transform → Digital HPF -----
z = poly(0, "z");
Hz = horner(H_hp_s, (2/T) * ((z-1)/(z+1)));

// ----- Frequency Response -----
Nf = 512;
HW = frmag(Hz, Nf);
w = 0 : %pi/(Nf-1) : %pi;

// ----- Plot -----
figure;
plot(w/%pi, abs(HW));
xlabel("Normalized Frequency (×π rad/sample)");
ylabel("Magnitude");
title("Chebyshev Type-I IIR High-Pass Filter Frequency Response");
xgrid();
```

## OUTPUT (LPF) : 

<img width="815" height="405" alt="3L" src="https://github.com/user-attachments/assets/1c5347a1-1dea-4cfd-abf8-c1c44ef9c350" />

## OUTPUT (HPF) : 

<img width="803" height="398" alt="3H" src="https://github.com/user-attachments/assets/99cf83d6-054c-4fba-b0f8-9ffbb0a87239" />


## RESULT: 

IIR Chebyshev filter(LPF & HPF) using SCILAB was designed successfully.
