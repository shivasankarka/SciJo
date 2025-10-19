"""
Collection of physical constants and conversion factors.

Author: Shivasankar K.A
Version: 0.1.0
Date: July 2025

Most constants are in SI units wherever applicable.
Based on SciPy constants module with 2022 CODATA values.
"""

import math

from numojo.core.datatypes import f64

# =======================
# MATHEMATICAL CONSTANTS
# =======================

alias pi: Scalar[f64] = 3.141592653589793  # math.pi
alias golden: Scalar[f64] = 1.618033988749895  # (1 + sqrt(5)) / 2
alias golden_ratio: Scalar[f64] = 1.618033988749895  # (1 + sqrt(5)) / 2

# =======================
# SI PREFIXES
# =======================

# Large prefixes
alias quetta: Scalar[f64] = 1e30
alias ronna: Scalar[f64] = 1e27
alias yotta: Scalar[f64] = 1e24
alias zetta: Scalar[f64] = 1e21
alias exa: Scalar[f64] = 1e18
alias peta: Scalar[f64] = 1e15
alias tera: Scalar[f64] = 1e12
alias giga: Scalar[f64] = 1e9
alias mega: Scalar[f64] = 1e6
alias kilo: Scalar[f64] = 1e3
alias hecto: Scalar[f64] = 1e2
alias deka: Scalar[f64] = 1e1

# Small prefixes
alias deci: Scalar[f64] = 1e-1
alias centi: Scalar[f64] = 1e-2
alias milli: Scalar[f64] = 1e-3
alias micro: Scalar[f64] = 1e-6
alias nano: Scalar[f64] = 1e-9
alias pico: Scalar[f64] = 1e-12
alias femto: Scalar[f64] = 1e-15
alias atto: Scalar[f64] = 1e-18
alias zepto: Scalar[f64] = 1e-21
alias yocto: Scalar[f64] = 1e-24
alias ronto: Scalar[f64] = 1e-27
alias quecto: Scalar[f64] = 1e-30

# =======================
# BINARY PREFIXES
# =======================

alias kibi: Scalar[f64] = 1024.0  # 2^10
alias mebi: Scalar[f64] = 1048576.0  # 2^20
alias gibi: Scalar[f64] = 1073741824.0  # 2^30
alias tebi: Scalar[f64] = 1099511627776.0  # 2^40
alias pebi: Scalar[f64] = 1125899906842624.0  # 2^50
alias exbi: Scalar[f64] = 1152921504606846976.0  # 2^60
alias zebi: Scalar[f64] = 1180591620717411303424.0  # 2^70
alias yobi: Scalar[f64] = 1208925819614629174706176.0  # 2^80

# =======================
# PHYSICAL CONSTANTS (2022 CODATA VALUES)
# =======================

# Speed of light
alias c: Scalar[f64] = 299792458.0  # speed of light in vacuum [m s^-1]
alias speed_of_light: Scalar[
    f64
] = 299792458.0  # speed of light in vacuum [m s^-1]

# Magnetic and electric constants
alias mu_0: Scalar[
    f64
] = 1.25663706127e-06  # vacuum magnetic permeability [N A^-2]
alias epsilon_0: Scalar[
    f64
] = 8.8541878188e-12  # vacuum electric permittivity [F m^-1]

# Planck constants
alias h: Scalar[f64] = 6.62607015e-34  # Planck constant [J Hz^-1]
alias Planck: Scalar[f64] = 6.62607015e-34  # Planck constant [J Hz^-1]
alias hbar: Scalar[
    f64
] = 1.0545718176461565e-34  # reduced Planck constant [J s]

# Gravitational constants
alias G: Scalar[
    f64
] = 6.6743e-11  # Newtonian constant of gravitation [m^3 kg^-1 s^-2]
alias gravitational_constant: Scalar[
    f64
] = 6.6743e-11  # Newtonian constant of gravitation [m^3 kg^-1 s^-2]
alias g: Scalar[f64] = 9.80665  # standard acceleration of gravity [m s^-2]

# Elementary charge
alias e: Scalar[f64] = 1.602176634e-19  # elementary charge [C]
alias elementary_charge: Scalar[f64] = 1.602176634e-19  # elementary charge [C]

# Gas constant
alias R: Scalar[f64] = 8.31446261815324  # molar gas constant [J mol^-1 K^-1]
alias gas_constant: Scalar[
    f64
] = 8.31446261815324  # molar gas constant [J mol^-1 K^-1]

# Fine structure constant
alias alpha: Scalar[
    f64
] = 0.0072973525643  # fine-structure constant [dimensionless]
alias fine_structure: Scalar[
    f64
] = 0.0072973525643  # fine-structure constant [dimensionless]

# Avogadro constant
alias N_A: Scalar[f64] = 6.02214076e23  # Avogadro constant [mol^-1]
alias Avogadro: Scalar[f64] = 6.02214076e23  # Avogadro constant [mol^-1]

# Boltzmann constant
alias k: Scalar[f64] = 1.380649e-23  # Boltzmann constant [J K^-1]
alias Boltzmann: Scalar[f64] = 1.380649e-23  # Boltzmann constant [J K^-1]

# Stefan-Boltzmann constant
alias sigma: Scalar[
    f64
] = 5.6703744191844314e-08  # Stefan-Boltzmann constant [W m^-2 K^-4]
alias Stefan_Boltzmann: Scalar[
    f64
] = 5.6703744191844314e-08  # Stefan-Boltzmann constant [W m^-2 K^-4]

# Wien displacement law constant
alias Wien: Scalar[
    f64
] = 0.0028977719551851727  # Wien wavelength displacement law constant [m K]

# Rydberg constant
alias Rydberg: Scalar[f64] = 10973731.568157  # Rydberg constant [m^-1]

# =======================
# MASS IN KG
# =======================

# Base mass units
alias gram: Scalar[f64] = 1e-3  # gram [kg]
alias metric_ton: Scalar[f64] = 1e3  # metric ton [kg]

# Imperial/US mass units
alias grain: Scalar[f64] = 64.79891e-6  # grain [kg]
alias lb: Scalar[f64] = 0.45359237  # pound (avoirdupois) [kg]
alias pound: Scalar[f64] = 0.45359237  # pound (avoirdupois) [kg]
alias oz: Scalar[f64] = 0.028349523125  # ounce [kg]:Scalar[f64] = pound/16
alias ounce: Scalar[f64] = 0.028349523125  # ounce [kg]:Scalar[f64] = pound/16
alias stone: Scalar[f64] = 6.35029318  # stone [kg]:Scalar[f64] = 14*pound
alias long_ton: Scalar[
    f64
] = 1016.0469088  # long ton [kg]:Scalar[f64] = 2240*pound
alias short_ton: Scalar[
    f64
] = 907.18474  # short ton [kg]:Scalar[f64] = 2000*pound

# Specialized mass units
alias troy_ounce: Scalar[
    f64
] = 0.0311034768  # troy ounce [kg]:Scalar[f64] = 480*grain
alias troy_pound: Scalar[
    f64
] = 0.3732417216  # troy pound [kg]:Scalar[f64] = 12*troy_ounce
alias carat: Scalar[f64] = 0.0002  # carat [kg]:Scalar[f64] = 200e-6
alias blob: Scalar[
    f64
] = 175.126835246  # blob [kg]:Scalar[f64] = pound*g/0.0254
alias slinch: Scalar[
    f64
] = 175.126835246  # slinch [kg]:Scalar[f64] = pound*g/0.0254
alias slug: Scalar[f64] = 14.593902937  # slug [kg]:Scalar[f64] = blob/12

# Particle masses (2022 CODATA values)
alias m_e: Scalar[f64] = 9.1093837139e-31  # electron mass [kg]
alias electron_mass: Scalar[f64] = 9.1093837139e-31  # electron mass [kg]
alias m_p: Scalar[f64] = 1.67262192595e-27  # proton mass [kg]
alias proton_mass: Scalar[f64] = 1.67262192595e-27  # proton mass [kg]
alias m_n: Scalar[f64] = 1.67492750056e-27  # neutron mass [kg]
alias neutron_mass: Scalar[f64] = 1.67492750056e-27  # neutron mass [kg]
alias m_u: Scalar[f64] = 1.66053906892e-27  # atomic mass constant [kg]
alias u: Scalar[f64] = 1.66053906892e-27  # atomic mass constant [kg]
alias atomic_mass: Scalar[f64] = 1.66053906892e-27  # atomic mass constant [kg]

# =======================
# ANGLE IN RADIANS
# =======================

alias degree: Scalar[
    f64
] = 0.017453292519943295  # degree [rad]:Scalar[f64] = pi/180
alias arcmin: Scalar[
    f64
] = 0.00029088820866572158  # arcminute [rad]:Scalar[f64] = degree/60
alias arcminute: Scalar[
    f64
] = 0.00029088820866572158  # arcminute [rad]:Scalar[f64] = degree/60
alias arcsec: Scalar[
    f64
] = 4.8481368110953599e-06  # arcsecond [rad]:Scalar[f64] = arcmin/60
alias arcsecond: Scalar[
    f64
] = 4.8481368110953599e-06  # arcsecond [rad]:Scalar[f64] = arcmin/60

# =======================
# TIME IN SECONDS
# =======================

alias minute: Scalar[f64] = 60.0  # minute [s]
alias hour: Scalar[f64] = 3600.0  # hour [s]:Scalar[f64] = 60*minute
alias day: Scalar[f64] = 86400.0  # day [s]:Scalar[f64] = 24*hour
alias week: Scalar[f64] = 604800.0  # week [s]:Scalar[f64] = 7*day
alias year: Scalar[f64] = 31536000.0  # year [s]:Scalar[f64] = 365*day
alias Julian_year: Scalar[
    f64
] = 31557600.0  # Julian year [s]:Scalar[f64] = 365.25*day

# =======================
# LENGTH IN METERS
# =======================

# Basic length units
alias inch: Scalar[f64] = 0.0254  # inch [m]
alias foot: Scalar[f64] = 0.3048  # foot [m]:Scalar[f64] = 12*inch
alias yard: Scalar[f64] = 0.9144  # yard [m]:Scalar[f64] = 3*foot
alias mile: Scalar[f64] = 1609.344  # mile [m]:Scalar[f64] = 1760*yard
alias mil: Scalar[f64] = 2.54e-05  # mil [m]:Scalar[f64] = inch/1000
alias pt: Scalar[
    f64
] = 0.00035277777777777776  # point [m]:Scalar[f64] = inch/72
alias point: Scalar[
    f64
] = 0.00035277777777777776  # point [m]:Scalar[f64] = inch/72

# Survey units
alias survey_foot: Scalar[
    f64
] = 0.30480060960121924  # survey foot [m]:Scalar[f64] = 1200.0/3937
alias survey_mile: Scalar[
    f64
] = 1609.3472186944375  # survey mile [m]:Scalar[f64] = 5280*survey_foot

# Maritime and scientific units
alias nautical_mile: Scalar[f64] = 1852.0  # nautical mile [m]
alias fermi: Scalar[f64] = 1e-15  # fermi [m]
alias angstrom: Scalar[f64] = 1e-10  # angstrom [m]
alias micron: Scalar[f64] = 1e-06  # micron [m]

# Astronomical units
alias au: Scalar[f64] = 149597870700.0  # astronomical unit [m]
alias astronomical_unit: Scalar[f64] = 149597870700.0  # astronomical unit [m]
alias light_year: Scalar[
    f64
] = 9460730472580800.0  # light year [m]:Scalar[f64] = Julian_year*c
alias parsec: Scalar[
    f64
] = 3.0856775814913673e16  # parsec [m]:Scalar[f64] = au/arcsec

# =======================
# PRESSURE IN PASCALS
# =======================

alias atm: Scalar[f64] = 101325.0  # standard atmosphere [Pa] (2022 CODATA)
alias atmosphere: Scalar[f64] = 101325.0  # standard atmosphere [Pa]
alias bar: Scalar[f64] = 100000.0  # bar [Pa]:Scalar[f64] = 1e5
alias torr: Scalar[f64] = 133.32236842105263  # torr [Pa]:Scalar[f64] = atm/760
alias mmHg: Scalar[f64] = 133.32236842105263  # mmHg [Pa]:Scalar[f64] = atm/760
alias psi: Scalar[
    f64
] = 6894.757293168361  # psi [Pa]:Scalar[f64] = pound*g/(inch*inch)

# =======================
# AREA IN SQUARE METERS
# =======================

alias hectare: Scalar[f64] = 10000.0  # hectare [m^2]:Scalar[f64] = 1e4
alias acre: Scalar[f64] = 4046.8564224  # acre [m^2]:Scalar[f64] = 43560*foot^2

# =======================
# VOLUME IN CUBIC METERS
# =======================

# Metric volume
alias liter: Scalar[f64] = 0.001  # liter [m^3]:Scalar[f64] = 1e-3
alias litre: Scalar[f64] = 0.001  # litre [m^3]:Scalar[f64] = 1e-3

# US volume
alias gallon: Scalar[
    f64
] = 0.003785411784  # gallon (US) [m^3]:Scalar[f64] = 231*inch^3
alias gallon_US: Scalar[
    f64
] = 0.003785411784  # gallon (US) [m^3]:Scalar[f64] = 231*inch^3
alias fluid_ounce: Scalar[
    f64
] = 2.95735295625e-05  # fluid ounce (US) [m^3]:Scalar[f64] = gallon_US/128
alias fluid_ounce_US: Scalar[
    f64
] = 2.95735295625e-05  # fluid ounce (US) [m^3]:Scalar[f64] = gallon_US/128
alias bbl: Scalar[
    f64
] = 0.158987294928  # barrel [m^3]:Scalar[f64] = 42*gallon_US
alias barrel: Scalar[
    f64
] = 0.158987294928  # barrel [m^3]:Scalar[f64] = 42*gallon_US

# UK volume
alias gallon_imp: Scalar[
    f64
] = 0.00454609  # gallon (UK) [m^3]:Scalar[f64] = 4.54609e-3
alias fluid_ounce_imp: Scalar[
    f64
] = 2.8413062499999997e-05  # fluid ounce (UK) [m^3]:Scalar[f64] = gallon_imp/160

# =======================
# SPEED IN METERS PER SECOND
# =======================

alias kmh: Scalar[
    f64
] = 0.2777777777777778  # km/h [m s^-1]:Scalar[f64] = 1e3/hour
alias mph: Scalar[f64] = 0.44704  # mph [m s^-1]:Scalar[f64] = mile/hour
alias mach: Scalar[f64] = 340.5  # Mach [m s^-1] (approx at 15°C, 1 atm)
alias speed_of_sound: Scalar[
    f64
] = 340.5  # speed of sound [m s^-1] (approx at 15°C, 1 atm)
alias knot: Scalar[
    f64
] = 0.5144444444444445  # knot [m s^-1]:Scalar[f64] = nautical_mile/hour

# =======================
# TEMPERATURE IN KELVIN
# =======================

alias zero_Celsius: Scalar[f64] = 273.15  # zero Celsius [K]
alias degree_Fahrenheit: Scalar[
    f64
] = 0.5555555555555556  # degree Fahrenheit [K] (for differences only):Scalar[f64] = 1/1.8

# =======================
# ENERGY IN JOULES
# =======================

# Basic energy units
alias eV: Scalar[
    f64
] = 1.602176634e-19  # electron volt [J]:Scalar[f64] = elementary_charge
alias electron_volt: Scalar[
    f64
] = 1.602176634e-19  # electron volt [J]:Scalar[f64] = elementary_charge

# Thermal energy units
alias calorie: Scalar[f64] = 4.184  # calorie (thermochemical) [J]
alias calorie_th: Scalar[f64] = 4.184  # calorie (thermochemical) [J]
alias calorie_IT: Scalar[
    f64
] = 4.1868  # calorie (International Steam Table) [J]
alias erg: Scalar[f64] = 1e-07  # erg [J]:Scalar[f64] = 1e-7

# British thermal units
alias Btu_th: Scalar[
    f64
] = 1054.3502644888888  # BTU (thermochemical) [J]:Scalar[f64] = pound*degree_Fahrenheit*calorie_th/gram
alias Btu: Scalar[
    f64
] = 1055.05585262  # BTU (International Steam Table) [J]:Scalar[f64] = pound*degree_Fahrenheit*calorie_IT/gram
alias Btu_IT: Scalar[f64] = 1055.05585262  # BTU (International Steam Table) [J]

# Explosive energy
alias ton_TNT: Scalar[
    f64
] = 4184000000.0  # ton of TNT [J]:Scalar[f64] = 1e9*calorie_th

# =======================
# POWER IN WATTS
# =======================

alias hp: Scalar[
    f64
] = 745.6998715822702  # horsepower [W]:Scalar[f64] = 550*foot*pound*g
alias horsepower: Scalar[
    f64
] = 745.6998715822702  # horsepower [W]:Scalar[f64] = 550*foot*pound*g

# =======================
# FORCE IN NEWTONS
# =======================

alias dyn: Scalar[f64] = 1e-05  # dyne [N]:Scalar[f64] = 1e-5
alias dyne: Scalar[f64] = 1e-05  # dyne [N]:Scalar[f64] = 1e-5
alias lbf: Scalar[
    f64
] = 4.4482216152605  # pound force [N]:Scalar[f64] = pound*g
alias pound_force: Scalar[
    f64
] = 4.4482216152605  # pound force [N]:Scalar[f64] = pound*g
alias kgf: Scalar[f64] = 9.80665  # kilogram force [N]:Scalar[f64] = g
alias kilogram_force: Scalar[
    f64
] = 9.80665  # kilogram force [N]:Scalar[f64] = g
