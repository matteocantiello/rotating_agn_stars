# MESA: Stellar Evolution in AGN Disks, including Rotation
Extensions to MESA (Modules for Experiments in Stellar Astrophysics) for modeling the evolution of stars embedded in Active Galactic Nucleus (AGN) disks.

## Overview

Stars embedded in AGN disks experience dramatically different conditions than isolated stars: high ambient densities, modified atmospheric boundary conditions, and continuous accretion of mass and angular momentum from the surrounding disk. This code extends MESA to self-consistently model these effects.

## Features

### Atmospheric Boundary Conditions
- Custom surface pressure and temperature boundary conditions appropriate for the AGN disk environment
- Gradual blending from isolated to embedded conditions to aid numerical stability

### Mass Accretion and Loss
- Bondi-like accretion from the AGN disk with radiative feedback
- Super-Eddington mass loss prescriptions
- Support for multiple accretion rate prescriptions (vertical rarefaction, shear, Hill, tidal barrier)

### Angular Momentum Evolution
Self-consistent angular momentum accretion following [Jermyn et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...914..105J), including:

- **Specific angular momentum of accreted material**: Computed as the minimum of the disk Keplerian value at the accretion radius and the stellar breakup limit:
  ```
  j_acc = min(R_acc² × Ω_AGN, √(GM★R★))
  ```
  where `R_acc = min(R_Hill, R_Bondi)`.

- **Angular momentum loss via mass loss**: Material ejected from the stellar surface carries away angular momentum at the local specific angular momentum:
  ```
  j_surf = Ω_surf × R★²
  ```

- **Net torque applied to the star**:
  ```
  dJ/dt = gain × j_acc - loss × j_surf
  ```

- **Rotational mixing**: Angular momentum is deposited throughout the star and MESA's angular momentum transport routines handle the internal redistribution, enabling self-consistent tracking of rotational mixing.

### Relevant Length Scales
- **Bondi radius**: `R_Bondi = 2GM★/c_s²`
- **Hill radius**: `R_Hill = (R_Bondi × h² / 12)^(1/3)` where `h` is the disk scale height
- **Accretion radius**: `R_acc = min(R_Hill, R_Bondi)`

## Requirements

- MESA r24.08.1

## Credits

Built on contributions from:
- Adam Jermyn
- Alex Dittmann  
- Ebraheem Farag
- Matteo Cantiello

## References

- Jermyn, A. S., Dittmann, A. J., Cantiello, M., & Perna, R. (2021). "Stellar Evolution in the Disks of Active Galactic Nuclei Produces Rapidly Rotating Massive Stars." *ApJ*, 914, 105. [ADS](https://ui.adsabs.harvard.edu/abs/2021ApJ...914..105J)

---

![PGPlot](rotating.png)

