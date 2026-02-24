"""
Collection of physical constants and conversion factors.

Author: Shivasankar K.A
Version: 0.1.0
Date: July 2025

Most constants are in SI units wherever applicable.
Based on SciPy constants module with 2022 CODATA values.
"""

import math

from numojo.core import f64

# =======================
# MATHEMATICAL CONSTANTS
# =======================

comptime pi: Scalar[f64] = 3.141592653589793  # math.pi
comptime golden: Scalar[f64] = 1.618033988749895  # (1 + sqrt(5)) / 2
comptime golden_ratio: Scalar[f64] = 1.618033988749895  # (1 + sqrt(5)) / 2

# =======================
# SI PREFIXES
# =======================

# Large prefixes
comptime quetta: Scalar[f64] = 1e30
comptime ronna: Scalar[f64] = 1e27
comptime yotta: Scalar[f64] = 1e24
comptime zetta: Scalar[f64] = 1e21
comptime exa: Scalar[f64] = 1e18
comptime peta: Scalar[f64] = 1e15
comptime tera: Scalar[f64] = 1e12
comptime giga: Scalar[f64] = 1e9
comptime mega: Scalar[f64] = 1e6
comptime kilo: Scalar[f64] = 1e3
comptime hecto: Scalar[f64] = 1e2
comptime deka: Scalar[f64] = 1e1

# Small prefixes
comptime deci: Scalar[f64] = 1e-1
comptime centi: Scalar[f64] = 1e-2
comptime milli: Scalar[f64] = 1e-3
comptime micro: Scalar[f64] = 1e-6
comptime nano: Scalar[f64] = 1e-9
comptime pico: Scalar[f64] = 1e-12
comptime femto: Scalar[f64] = 1e-15
comptime atto: Scalar[f64] = 1e-18
comptime zepto: Scalar[f64] = 1e-21
comptime yocto: Scalar[f64] = 1e-24
comptime ronto: Scalar[f64] = 1e-27
comptime quecto: Scalar[f64] = 1e-30

# =======================
# BINARY PREFIXES
# =======================

comptime kibi: Scalar[f64] = 1024.0  # 2^10
comptime mebi: Scalar[f64] = 1048576.0  # 2^20
comptime gibi: Scalar[f64] = 1073741824.0  # 2^30
comptime tebi: Scalar[f64] = 1099511627776.0  # 2^40
comptime pebi: Scalar[f64] = 1125899906842624.0  # 2^50
comptime exbi: Scalar[f64] = 1152921504606846976.0  # 2^60
comptime zebi: Scalar[f64] = 1180591620717411303424.0  # 2^70
comptime yobi: Scalar[f64] = 1208925819614629174706176.0  # 2^80

# =======================
# PHYSICAL CONSTANTS (2022 CODATA VALUES)
# =======================

# Speed of light
comptime c: Scalar[f64] = 299792458.0  # speed of light in vacuum [m s^-1]
comptime speed_of_light: Scalar[
    f64
] = 299792458.0  # speed of light in vacuum [m s^-1]

# Magnetic and electric constants
comptime mu_0: Scalar[
    f64
] = 1.25663706127e-06  # vacuum magnetic permeability [N A^-2]
comptime epsilon_0: Scalar[
    f64
] = 8.8541878188e-12  # vacuum electric permittivity [F m^-1]

# Planck constants
comptime h: Scalar[f64] = 6.62607015e-34  # Planck constant [J Hz^-1]
comptime Planck: Scalar[f64] = 6.62607015e-34  # Planck constant [J Hz^-1]
comptime hbar: Scalar[
    f64
] = 1.0545718176461565e-34  # reduced Planck constant [J s]

# Gravitational constants
comptime G: Scalar[
    f64
] = 6.6743e-11  # Newtonian constant of gravitation [m^3 kg^-1 s^-2]
comptime gravitational_constant: Scalar[
    f64
] = 6.6743e-11  # Newtonian constant of gravitation [m^3 kg^-1 s^-2]
comptime g: Scalar[f64] = 9.80665  # standard acceleration of gravity [m s^-2]

# Elementary charge
comptime e: Scalar[f64] = 1.602176634e-19  # elementary charge [C]
comptime elementary_charge: Scalar[
    f64
] = 1.602176634e-19  # elementary charge [C]

# Gas constant
comptime R: Scalar[f64] = 8.31446261815324  # molar gas constant [J mol^-1 K^-1]
comptime gas_constant: Scalar[
    f64
] = 8.31446261815324  # molar gas constant [J mol^-1 K^-1]

# Fine structure constant
comptime alpha: Scalar[
    f64
] = 0.0072973525643  # fine-structure constant [dimensionless]
comptime fine_structure: Scalar[
    f64
] = 0.0072973525643  # fine-structure constant [dimensionless]

# Avogadro constant
comptime N_A: Scalar[f64] = 6.02214076e23  # Avogadro constant [mol^-1]
comptime Avogadro: Scalar[f64] = 6.02214076e23  # Avogadro constant [mol^-1]

# Boltzmann constant
comptime k: Scalar[f64] = 1.380649e-23  # Boltzmann constant [J K^-1]
comptime Boltzmann: Scalar[f64] = 1.380649e-23  # Boltzmann constant [J K^-1]

# Stefan-Boltzmann constant
comptime sigma: Scalar[
    f64
] = 5.6703744191844314e-08  # Stefan-Boltzmann constant [W m^-2 K^-4]
comptime Stefan_Boltzmann: Scalar[
    f64
] = 5.6703744191844314e-08  # Stefan-Boltzmann constant [W m^-2 K^-4]

# Wien displacement law constant
comptime Wien: Scalar[
    f64
] = 0.0028977719551851727  # Wien wavelength displacement law constant [m K]

# Rydberg constant
comptime Rydberg: Scalar[f64] = 10973731.568157  # Rydberg constant [m^-1]

# =======================
# MASS IN KG
# =======================

# Base mass units
comptime gram: Scalar[f64] = 1e-3  # gram [kg]
comptime metric_ton: Scalar[f64] = 1e3  # metric ton [kg]

# Imperial/US mass units
comptime grain: Scalar[f64] = 64.79891e-6  # grain [kg]
comptime lb: Scalar[f64] = 0.45359237  # pound (avoirdupois) [kg]
comptime pound: Scalar[f64] = 0.45359237  # pound (avoirdupois) [kg]
comptime oz: Scalar[f64] = 0.028349523125  # ounce [kg]:Scalar[f64] = pound/16
comptime ounce: Scalar[
    f64
] = 0.028349523125  # ounce [kg]:Scalar[f64] = pound/16
comptime stone: Scalar[f64] = 6.35029318  # stone [kg]:Scalar[f64] = 14*pound
comptime long_ton: Scalar[
    f64
] = 1016.0469088  # long ton [kg]:Scalar[f64] = 2240*pound
comptime short_ton: Scalar[
    f64
] = 907.18474  # short ton [kg]:Scalar[f64] = 2000*pound

# Specialized mass units
comptime troy_ounce: Scalar[
    f64
] = 0.0311034768  # troy ounce [kg]:Scalar[f64] = 480*grain
comptime troy_pound: Scalar[
    f64
] = 0.3732417216  # troy pound [kg]:Scalar[f64] = 12*troy_ounce
comptime carat: Scalar[f64] = 0.0002  # carat [kg]:Scalar[f64] = 200e-6
comptime blob: Scalar[
    f64
] = 175.126835246  # blob [kg]:Scalar[f64] = pound*g/0.0254
comptime slinch: Scalar[
    f64
] = 175.126835246  # slinch [kg]:Scalar[f64] = pound*g/0.0254
comptime slug: Scalar[f64] = 14.593902937  # slug [kg]:Scalar[f64] = blob/12

# Particle masses (2022 CODATA values)
comptime m_e: Scalar[f64] = 9.1093837139e-31  # electron mass [kg]
comptime electron_mass: Scalar[f64] = 9.1093837139e-31  # electron mass [kg]
comptime m_p: Scalar[f64] = 1.67262192595e-27  # proton mass [kg]
comptime proton_mass: Scalar[f64] = 1.67262192595e-27  # proton mass [kg]
comptime m_n: Scalar[f64] = 1.67492750056e-27  # neutron mass [kg]
comptime neutron_mass: Scalar[f64] = 1.67492750056e-27  # neutron mass [kg]
comptime m_u: Scalar[f64] = 1.66053906892e-27  # atomic mass constant [kg]
comptime u: Scalar[f64] = 1.66053906892e-27  # atomic mass constant [kg]
comptime atomic_mass: Scalar[
    f64
] = 1.66053906892e-27  # atomic mass constant [kg]

# =======================
# ANGLE IN RADIANS
# =======================

comptime degree: Scalar[
    f64
] = 0.017453292519943295  # degree [rad]:Scalar[f64] = pi/180
comptime arcmin: Scalar[
    f64
] = 0.00029088820866572158  # arcminute [rad]:Scalar[f64] = degree/60
comptime arcminute: Scalar[
    f64
] = 0.00029088820866572158  # arcminute [rad]:Scalar[f64] = degree/60
comptime arcsec: Scalar[
    f64
] = 4.8481368110953599e-06  # arcsecond [rad]:Scalar[f64] = arcmin/60
comptime arcsecond: Scalar[
    f64
] = 4.8481368110953599e-06  # arcsecond [rad]:Scalar[f64] = arcmin/60

# =======================
# TIME IN SECONDS
# =======================

comptime minute: Scalar[f64] = 60.0  # minute [s]
comptime hour: Scalar[f64] = 3600.0  # hour [s]:Scalar[f64] = 60*minute
comptime day: Scalar[f64] = 86400.0  # day [s]:Scalar[f64] = 24*hour
comptime week: Scalar[f64] = 604800.0  # week [s]:Scalar[f64] = 7*day
comptime year: Scalar[f64] = 31536000.0  # year [s]:Scalar[f64] = 365*day
comptime Julian_year: Scalar[
    f64
] = 31557600.0  # Julian year [s]:Scalar[f64] = 365.25*day

# =======================
# LENGTH IN METERS
# =======================

# Basic length units
comptime inch: Scalar[f64] = 0.0254  # inch [m]
comptime foot: Scalar[f64] = 0.3048  # foot [m]:Scalar[f64] = 12*inch
comptime yard: Scalar[f64] = 0.9144  # yard [m]:Scalar[f64] = 3*foot
comptime mile: Scalar[f64] = 1609.344  # mile [m]:Scalar[f64] = 1760*yard
comptime mil: Scalar[f64] = 2.54e-05  # mil [m]:Scalar[f64] = inch/1000
comptime pt: Scalar[
    f64
] = 0.00035277777777777776  # point [m]:Scalar[f64] = inch/72
comptime point: Scalar[
    f64
] = 0.00035277777777777776  # point [m]:Scalar[f64] = inch/72

# Survey units
comptime survey_foot: Scalar[
    f64
] = 0.30480060960121924  # survey foot [m]:Scalar[f64] = 1200.0/3937
comptime survey_mile: Scalar[
    f64
] = 1609.3472186944375  # survey mile [m]:Scalar[f64] = 5280*survey_foot

# Maritime and scientific units
comptime nautical_mile: Scalar[f64] = 1852.0  # nautical mile [m]
comptime fermi: Scalar[f64] = 1e-15  # fermi [m]
comptime angstrom: Scalar[f64] = 1e-10  # angstrom [m]
comptime micron: Scalar[f64] = 1e-06  # micron [m]

# Astronomical units
comptime au: Scalar[f64] = 149597870700.0  # astronomical unit [m]
comptime astronomical_unit: Scalar[
    f64
] = 149597870700.0  # astronomical unit [m]
comptime light_year: Scalar[
    f64
] = 9460730472580800.0  # light year [m]:Scalar[f64] = Julian_year*c
comptime parsec: Scalar[
    f64
] = 3.0856775814913673e16  # parsec [m]:Scalar[f64] = au/arcsec

# =======================
# PRESSURE IN PASCALS
# =======================

comptime atm: Scalar[f64] = 101325.0  # standard atmosphere [Pa] (2022 CODATA)
comptime atmosphere: Scalar[f64] = 101325.0  # standard atmosphere [Pa]
comptime bar: Scalar[f64] = 100000.0  # bar [Pa]:Scalar[f64] = 1e5
comptime torr: Scalar[
    f64
] = 133.32236842105263  # torr [Pa]:Scalar[f64] = atm/760
comptime mmHg: Scalar[
    f64
] = 133.32236842105263  # mmHg [Pa]:Scalar[f64] = atm/760
comptime psi: Scalar[
    f64
] = 6894.757293168361  # psi [Pa]:Scalar[f64] = pound*g/(inch*inch)

# =======================
# AREA IN SQUARE METERS
# =======================

comptime hectare: Scalar[f64] = 10000.0  # hectare [m^2]:Scalar[f64] = 1e4
comptime acre: Scalar[
    f64
] = 4046.8564224  # acre [m^2]:Scalar[f64] = 43560*foot^2

# =======================
# VOLUME IN CUBIC METERS
# =======================

# Metric volume
comptime liter: Scalar[f64] = 0.001  # liter [m^3]:Scalar[f64] = 1e-3
comptime litre: Scalar[f64] = 0.001  # litre [m^3]:Scalar[f64] = 1e-3

# US volume
comptime gallon: Scalar[
    f64
] = 0.003785411784  # gallon (US) [m^3]:Scalar[f64] = 231*inch^3
comptime gallon_US: Scalar[
    f64
] = 0.003785411784  # gallon (US) [m^3]:Scalar[f64] = 231*inch^3
comptime fluid_ounce: Scalar[
    f64
] = 2.95735295625e-05  # fluid ounce (US) [m^3]:Scalar[f64] = gallon_US/128
comptime fluid_ounce_US: Scalar[
    f64
] = 2.95735295625e-05  # fluid ounce (US) [m^3]:Scalar[f64] = gallon_US/128
comptime bbl: Scalar[
    f64
] = 0.158987294928  # barrel [m^3]:Scalar[f64] = 42*gallon_US
comptime barrel: Scalar[
    f64
] = 0.158987294928  # barrel [m^3]:Scalar[f64] = 42*gallon_US

# UK volume
comptime gallon_imp: Scalar[
    f64
] = 0.00454609  # gallon (UK) [m^3]:Scalar[f64] = 4.54609e-3
comptime fluid_ounce_imp: Scalar[
    f64
] = 2.8413062499999997e-05  # fluid ounce (UK) [m^3]:Scalar[f64] = gallon_imp/160

# =======================
# SPEED IN METERS PER SECOND
# =======================

comptime kmh: Scalar[
    f64
] = 0.2777777777777778  # km/h [m s^-1]:Scalar[f64] = 1e3/hour
comptime mph: Scalar[f64] = 0.44704  # mph [m s^-1]:Scalar[f64] = mile/hour
comptime mach: Scalar[f64] = 340.5  # Mach [m s^-1] (approx at 15°C, 1 atm)
comptime speed_of_sound: Scalar[
    f64
] = 340.5  # speed of sound [m s^-1] (approx at 15°C, 1 atm)
comptime knot: Scalar[
    f64
] = 0.5144444444444445  # knot [m s^-1]:Scalar[f64] = nautical_mile/hour

# =======================
# TEMPERATURE IN KELVIN
# =======================

comptime zero_Celsius: Scalar[f64] = 273.15  # zero Celsius [K]
comptime degree_Fahrenheit: Scalar[
    f64
] = 0.5555555555555556  # degree Fahrenheit [K] (for differences only):Scalar[f64] = 1/1.8

# =======================
# ENERGY IN JOULES
# =======================

# Basic energy units
comptime eV: Scalar[
    f64
] = 1.602176634e-19  # electron volt [J]:Scalar[f64] = elementary_charge
comptime electron_volt: Scalar[
    f64
] = 1.602176634e-19  # electron volt [J]:Scalar[f64] = elementary_charge

# Thermal energy units
comptime calorie: Scalar[f64] = 4.184  # calorie (thermochemical) [J]
comptime calorie_th: Scalar[f64] = 4.184  # calorie (thermochemical) [J]
comptime calorie_IT: Scalar[
    f64
] = 4.1868  # calorie (International Steam Table) [J]
comptime erg: Scalar[f64] = 1e-07  # erg [J]:Scalar[f64] = 1e-7

# British thermal units
comptime Btu_th: Scalar[
    f64
] = 1054.3502644888888  # BTU (thermochemical) [J]:Scalar[f64] = pound*degree_Fahrenheit*calorie_th/gram
comptime Btu: Scalar[
    f64
] = 1055.05585262  # BTU (International Steam Table) [J]:Scalar[f64] = pound*degree_Fahrenheit*calorie_IT/gram
comptime Btu_IT: Scalar[
    f64
] = 1055.05585262  # BTU (International Steam Table) [J]

# Explosive energy
comptime ton_TNT: Scalar[
    f64
] = 4184000000.0  # ton of TNT [J]:Scalar[f64] = 1e9*calorie_th

# =======================
# POWER IN WATTS
# =======================

comptime hp: Scalar[
    f64
] = 745.6998715822702  # horsepower [W]:Scalar[f64] = 550*foot*pound*g
comptime horsepower: Scalar[
    f64
] = 745.6998715822702  # horsepower [W]:Scalar[f64] = 550*foot*pound*g

# =======================
# FORCE IN NEWTONS
# =======================

comptime dyn: Scalar[f64] = 1e-05  # dyne [N]:Scalar[f64] = 1e-5
comptime dyne: Scalar[f64] = 1e-05  # dyne [N]:Scalar[f64] = 1e-5
comptime lbf: Scalar[
    f64
] = 4.4482216152605  # pound force [N]:Scalar[f64] = pound*g
comptime pound_force: Scalar[
    f64
] = 4.4482216152605  # pound force [N]:Scalar[f64] = pound*g
comptime kgf: Scalar[f64] = 9.80665  # kilogram force [N]:Scalar[f64] = g
comptime kilogram_force: Scalar[
    f64
] = 9.80665  # kilogram force [N]:Scalar[f64] = g
