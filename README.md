**English** | [日本語](README_JP.md)

# ALICE-Chemistry

Molecular dynamics and chemistry simulation library for [Project A.L.I.C.E.](https://github.com/anthropics/alice)

## Overview

`alice-chemistry` provides a pure Rust chemistry simulation toolkit covering molecular dynamics, reaction kinetics, thermodynamics, and the periodic table.

## Features

- **Periodic Table** — element data for H through Kr (atomic number, mass, electronegativity)
- **Force Fields** — Lennard-Jones and Coulomb potentials
- **Reaction Kinetics** — Arrhenius equation rate constant calculation
- **Chemical Equilibrium** — equilibrium constant and concentration computation
- **Molecular Structure** — bond energy and molecular geometry
- **Thermodynamics** — enthalpy, entropy, Gibbs free energy calculations
- **Stoichiometry** — balanced reaction coefficient handling
- **Physical Constants** — Boltzmann, Avogadro, gas constant, Coulomb constant, etc.

## Quick Start

```rust
use alice_chemistry::{lookup_element, BOLTZMANN, GAS_CONSTANT};

let hydrogen = lookup_element(1).unwrap();
assert_eq!(hydrogen.symbol, "H");
assert!((hydrogen.atomic_mass - 1.008).abs() < 0.001);

// Arrhenius rate constant: k = A * exp(-Ea / (R * T))
let k = 1e13_f64 * (-75000.0 / (GAS_CONSTANT * 298.15)).exp();
```

## Architecture

```
alice-chemistry
├── Element           — element data (atomic number, mass, electronegativity)
├── ELEMENTS          — periodic table data (H-Kr)
├── lennard_jones     — LJ potential & force calculation
├── coulomb           — electrostatic potential
├── arrhenius         — reaction rate constant
├── equilibrium       — chemical equilibrium solver
├── thermodynamics    — enthalpy, entropy, Gibbs energy
└── stoichiometry     — reaction balancing & molar ratios
```

## License

MIT OR Apache-2.0
