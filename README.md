# LI‑COR Methane Flux Shiny App

This repository contains a Shiny application to process and calculate methane (CH₄) fluxes from static floating chamber measurements connected to a LI‑COR portable greenhouse gas analyzer (LI‑7810).

## Overview

Objective: Compute methane fluxes (mg m⁻² s⁻¹) using closed‑loop chamber data.

Method: Adapted from Gerardo‑Nieto et al. (2017), using a 4.24 L chamber with 453 cm² surface area. Fluxes are derived from the slope of the linear regression (ΔC/Δt) of concentration vs. time.

Data: Time series of CH₄ concentration recorded at 1 Hz. First 30 s discarded. Only regressions with R² > 0.82 are accepted. All measurements in triplicate.

### Contents
```bash
flux_app/
├── app.R             # Main Shiny application
└── data/             # Example dataset
    ├── all_data.csv   # Raw gas concentration time series
    └── list_all.csv   # Metadata list of measurements
```
### Installation

Ensure you have R (≥ 4.0) installed. Then install required packages:
```bash
install.packages(c("shiny", "ggplot2", "DT", "here"))
```
### Usage

From the project root directory:

# Launch the Shiny app:
shiny::runApp("flux_app")

In the sidebar, click Next sample to iterate over each measurement in list_all.csv.

Adjust Seconds before/after and Zoom‑out padding for CH₄ time series.

View the concentration vs. time plot, regression outputs, and computed flux (mg m⁻² s⁻¹).

Click Save computed fluxes to store results in a table and export to your R global environment as df_flux_final.

# Method Details

## Chamber setup (user-configurable)

Default volume: 4.24 L (set voltotal in the app header)

Default surface area: 453 cm² (set surface_area in the app header)

## Sampling protocol

CH₄ concentration measured at 1 Hz.

It is recommended to discard the first 30 s (configurable via the seconds-before setting) to avoid deployment disturbances.

## Regression criteria

Linear fit of concentration vs. time over 1–3 min

# Flux calculation:

F = (ΔC/Δt) × (Vc / Ac) × Mg

F: flux (mg m⁻² s⁻¹)

ΔC/Δt: slope from regression (volume fraction per second)

Vc: chamber volume (L)

Ac: chamber area (m²)

Mg: molecular mass of CH₄ (16.04 g mol⁻¹)

# Citation

If you use this code, please cite the associated article:

Author(s) (2025). Title. Journal, vol, pages. DOI:XXXX
