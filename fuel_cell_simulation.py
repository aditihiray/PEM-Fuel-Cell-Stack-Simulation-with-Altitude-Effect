# Fuel Cell Stack Simulation Considering Altitude Effects
# Author: Aditi Atul Hiray
# Based on: IMEKO Conference Paper

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Physical Constants
# ----------------------------
R = 8.314       # J/(mol·K), Universal gas constant
F = 96485       # C/mol, Faraday's constant
T0 = 288.15     # Sea level temperature in K
P0 = 101.3      # Sea level pressure in kPa
g = 9.80665     # m/s²
M = 0.0289644   # kg/mol, molar mass of air
L = 0.0065      # K/m, lapse rate
Tref = 298.15   # Reference temperature in K
n = 2           # Electrons per reaction
E0 = 1.229      # Standard reversible potential at STP [V]
delta_S = -163  # Entropy change [J/mol·K]

# ----------------------------
# Fuel Cell Parameters
# ----------------------------
alpha = 0.5     # Electron transfer coefficient
i0 = 1e-4       # Exchange current density [A/cm²]
iL = 0.81       # Limiting current density [A/cm²]
kc = 1e-3       # Reaction rate factor
beta = 0.5      # Symmetry factor
tm = 0.0127     # Membrane thickness [cm]
ld = 0.03       # Diffusion layer thickness [cm]
sigma_d = 5     # Diffusion layer conductivity [S/cm]
phi = 0.95      # Relative humidity
i_internal = 0.002  # Internal loss current [A/cm²]
A = 50          # Active area [cm²]
n_cells = 10    # Cells in stack

# ----------------------------
# Altitude Models
# ----------------------------
def temperature_at_altitude(h):
    return T0 - L * h

def pressure_at_altitude(h):
    T = temperature_at_altitude(h)
    return P0 * (T / T0) ** (g * M / (R * L))

def oxygen_partial_pressure(h):
    return 0.21 * pressure_at_altitude(h)

# ----------------------------
# Membrane Water Content
# ----------------------------
def membrane_water_content(phi):
    return 0.043 + 17.81 * phi - 39.85 * phi**2 + 36 * phi**3

# ----------------------------
# Ionic Resistance
# ----------------------------
def ionic_resistance(T, i):
    lam = membrane_water_content(phi)
    num = 181.6 * (1 + 0.03 * i + 0.062 * i**2)
    denom = (T - 303) * lam - 0.634 - 3 * np.exp(4.18 * ((T - 303) / T))
    return tm * num / (303 * denom + 1e-6)  # add epsilon to avoid div by 0

# ----------------------------
# Voltage Loss Terms
# ----------------------------
def reversible_voltage(T, p_O2, p_H2, p_H2O):
    p_O2 = max(p_O2, 1e-3)
    p_H2 = max(p_H2, 1e-3)
    p_H2O = max(p_H2O, 1e-3)
    term1 = E0
    term2 = -delta_S * (T - Tref) / (n * F)
    term3 = (R * T) / (2 * n * F) * np.log(p_H2 * np.sqrt(p_O2) / p_H2O)
    return term1 + term2 + term3

def activation_loss(T, i):
    i = max(i, 1e-6)
    i0_T = kc * np.exp(-n * F * beta / (R * T))  # i0 ~ f(T)
    return (R * T) / (alpha * n * F) * np.log(i / i0_T)

def concentration_loss(T, i):
    i_eff = i + i_internal
    ratio = np.clip(i_eff / iL, 1e-6, 0.999)
    return -(R * Tref) / (n * F) * np.log(1 - ratio)

def ohmic_loss(T, i):
    r_ion = ionic_resistance(T, i)
    r_el = ld**2 / sigma_d
    return i * (r_ion + r_el)

# ----------------------------
# Total Fuel Cell Voltage
# ----------------------------
def fuel_cell_voltage(i, T, p_O2, p_H2=1.0, p_H2O=0.0313):
    E = reversible_voltage(T, p_O2, p_H2, p_H2O)
    V_act = activation_loss(T, i)
    V_conc = concentration_loss(T, i)
    V_ohm = ohmic_loss(T, i)
    V_cell = E - V_act - V_conc - V_ohm
    return max(V_cell, 0)

# ----------------------------
# Simulation Loop
# ----------------------------
altitudes = [0, 2000, 4000, 6000, 8000]  # meters
I = np.linspace(0.01, 1.0, 100)          # stack current in A
i_density = I / A                        # current density [A/cm²]

plt.figure(figsize=(9, 5))

for h in altitudes:
    T = temperature_at_altitude(h)
    p_O2 = oxygen_partial_pressure(h)
    V_stack = []
    for i in i_density:
        V_cell = fuel_cell_voltage(i, T, p_O2)
        V_stack.append(V_cell * n_cells)
    plt.plot(I, V_stack, label=f'{h} m, {round(T - 273.15, 1)} °C')

# ----------------------------
# Plot Settings
# ----------------------------
plt.title('Fuel Cell Stack Voltage vs. Current')
plt.xlabel('Stack Current (A)')
plt.ylabel('Stack Voltage (V)')
plt.legend(title='Altitude & Temperature')
plt.grid(True)
plt.tight_layout()
plt.show()
