# SIDARTHE Epidemic Model in Python

A numerical implementation of the SIDARTHE compartmental model for COVID-19, from Giordano et al., *"Modelling the COVID-19 epidemic and implementation of population-wide interventions in Italy"*, Nature Medicine 26, 855-860 (2020).

*Author:* [Dezső Szilárd, Speilmann Patrik]

## Introduction

SIDARTHE extends the classic SIR model by separating infected individuals according to whether they have been **detected** and how **severe** their symptoms are. 

The population is split into eight parts:

| Symbol | Meaning |
|--------|---------|
| S | Susceptible |
| I | Infected: asymptomatic, undetected |
| D | Diagnosed: asymptomatic, detected |
| A | Ailing: symptomatic, undetected |
| R | Recognized: symptomatic, detected |
| T | Threatened: life-threatening symptoms, detected |
| H | Healed |
| E | Extinct (deceased) |

The dynamics are eight coupled ODEs. Transmission comes from four infectious compartments (I, D, A, R) with rates α, β, γ, δ. Detection, symptom progression, recovery and mortality are governed by the remaining parameters (ε, θ, ζ, η, μ, ν, τ, λ, κ, ξ, ρ, σ).

The basic reproduction number is a sum over transmission pathways:

R₀ = α/r₁ + βε/(r₁r₂) + γζ/(r₁r₃) + δηε/(r₁r₂r₄) + δζθ/(r₁r₃r₄)

with r₁ = ε+ζ+λ, r₂ = η+ρ, r₃ = θ+μ+κ, r₄ = ν+ξ. The paper reports R₀ = 2.38 for the initial parameter set.

## Goals

1. Implement the model and verify R₀ ≈ 2.38 with the paper's parameters.
2. Simulate the baseline epidemic and validate it (population conservation, early growth rate vs. R₀).
3. Simulate interventions by changing parameters: social distancing (lower α, β, γ, δ), lockdown (stronger reductions), and increased testing (higher ε, θ).
4. Compare scenarios: peak size and timing, total deaths and recoveries.
5. (Stretch) Time-varying parameters, vaccination or partial immunity.

## Setup

Parameter values and intervention schedule follow the paper's Methods section.
Initial conditions: S(0) = 0.999, I(0) = 0.001, all others 0. 

