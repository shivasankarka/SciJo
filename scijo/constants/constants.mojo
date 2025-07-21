"""
Collection of physical constants and conversion factors.

Author: Shivasankar K.A
Version: 0.1.0
Date: July 2025

Most constants are in SI units, so you can do calculations like:
print "10 mile per minute is", 10*mile/minute, "m/s"
or print "10 mile per minute is", 10*mile/(minute*knot), "knots"

Based on SciPy constants module with 2022 CODATA values.
"""

import math

# =======================
# MATHEMATICAL CONSTANTS
# =======================

alias pi = 3.141592653589793  # math.pi equivalent
alias golden = 1.618033988749895  # (1 + sqrt(5)) / 2
alias golden_ratio = 1.618033988749895  # (1 + sqrt(5)) / 2

# =======================
# SI PREFIXES
# =======================

# Large prefixes
alias quetta = 1e30
alias ronna = 1e27
alias yotta = 1e24
alias zetta = 1e21
alias exa = 1e18
alias peta = 1e15
alias tera = 1e12
alias giga = 1e9
alias mega = 1e6
alias kilo = 1e3
alias hecto = 1e2
alias deka = 1e1

# Small prefixes
alias deci = 1e-1
alias centi = 1e-2
alias milli = 1e-3
alias micro = 1e-6
alias nano = 1e-9
alias pico = 1e-12
alias femto = 1e-15
alias atto = 1e-18
alias zepto = 1e-21
alias yocto = 1e-24
alias ronto = 1e-27
alias quecto = 1e-30

# =======================
# BINARY PREFIXES
# =======================

alias kibi = 1024.0  # 2^10
alias mebi = 1048576.0  # 2^20
alias gibi = 1073741824.0  # 2^30
alias tebi = 1099511627776.0  # 2^40
alias pebi = 1125899906842624.0  # 2^50
alias exbi = 1152921504606846976.0  # 2^60
alias zebi = 1180591620717411303424.0  # 2^70
alias yobi = 1208925819614629174706176.0  # 2^80

# =======================
# PHYSICAL CONSTANTS (2022 CODATA VALUES)
# =======================

# Speed of light
alias c = 299792458.0  # speed of light in vacuum [m s^-1]
alias speed_of_light = 299792458.0  # speed of light in vacuum [m s^-1]

# Magnetic and electric constants
alias mu_0 = 1.25663706127e-06  # vacuum magnetic permeability [N A^-2]
alias epsilon_0 = 8.8541878188e-12  # vacuum electric permittivity [F m^-1]

# Planck constants
alias h = 6.62607015e-34  # Planck constant [J Hz^-1]
alias Planck = 6.62607015e-34  # Planck constant [J Hz^-1]
alias hbar = 1.0545718176461565e-34  # reduced Planck constant [J s]

# Gravitational constants
alias G = 6.6743e-11  # Newtonian constant of gravitation [m^3 kg^-1 s^-2]
alias gravitational_constant = 6.6743e-11  # Newtonian constant of gravitation [m^3 kg^-1 s^-2]
alias g = 9.80665  # standard acceleration of gravity [m s^-2]

# Elementary charge
alias e = 1.602176634e-19  # elementary charge [C]
alias elementary_charge = 1.602176634e-19  # elementary charge [C]

# Gas constant
alias R = 8.31446261815324  # molar gas constant [J mol^-1 K^-1]
alias gas_constant = 8.31446261815324  # molar gas constant [J mol^-1 K^-1]

# Fine structure constant
alias alpha = 0.0072973525643  # fine-structure constant [dimensionless]
alias fine_structure = 0.0072973525643  # fine-structure constant [dimensionless]

# Avogadro constant
alias N_A = 6.02214076e23  # Avogadro constant [mol^-1]
alias Avogadro = 6.02214076e23  # Avogadro constant [mol^-1]

# Boltzmann constant
alias k = 1.380649e-23  # Boltzmann constant [J K^-1]
alias Boltzmann = 1.380649e-23  # Boltzmann constant [J K^-1]

# Stefan-Boltzmann constant
alias sigma = 5.6703744191844314e-08  # Stefan-Boltzmann constant [W m^-2 K^-4]
alias Stefan_Boltzmann = 5.6703744191844314e-08  # Stefan-Boltzmann constant [W m^-2 K^-4]

# Wien displacement law constant
alias Wien = 0.0028977719551851727  # Wien wavelength displacement law constant [m K]

# Rydberg constant
alias Rydberg = 10973731.568157  # Rydberg constant [m^-1]

# =======================
# MASS IN KG
# =======================

# Base mass units
alias gram = 1e-3  # gram [kg]
alias metric_ton = 1e3  # metric ton [kg]

# Imperial/US mass units
alias grain = 64.79891e-6  # grain [kg]
alias lb = 0.45359237  # pound (avoirdupois) [kg]
alias pound = 0.45359237  # pound (avoirdupois) [kg]
alias oz = 0.028349523125  # ounce [kg] = pound/16
alias ounce = 0.028349523125  # ounce [kg] = pound/16
alias stone = 6.35029318  # stone [kg] = 14*pound
alias long_ton = 1016.0469088  # long ton [kg] = 2240*pound
alias short_ton = 907.18474  # short ton [kg] = 2000*pound

# Specialized mass units
alias troy_ounce = 0.0311034768  # troy ounce [kg] = 480*grain
alias troy_pound = 0.3732417216  # troy pound [kg] = 12*troy_ounce
alias carat = 0.0002  # carat [kg] = 200e-6
alias blob = 175.126835246  # blob [kg] = pound*g/0.0254
alias slinch = 175.126835246  # slinch [kg] = pound*g/0.0254
alias slug = 14.593902937  # slug [kg] = blob/12

# Particle masses (2022 CODATA values)
alias m_e = 9.1093837139e-31  # electron mass [kg]
alias electron_mass = 9.1093837139e-31  # electron mass [kg]
alias m_p = 1.67262192595e-27  # proton mass [kg]
alias proton_mass = 1.67262192595e-27  # proton mass [kg]
alias m_n = 1.67492750056e-27  # neutron mass [kg]
alias neutron_mass = 1.67492750056e-27  # neutron mass [kg]
alias m_u = 1.66053906892e-27  # atomic mass constant [kg]
alias u = 1.66053906892e-27  # atomic mass constant [kg]
alias atomic_mass = 1.66053906892e-27  # atomic mass constant [kg]

# =======================
# ANGLE IN RADIANS
# =======================

alias degree = 0.017453292519943295  # degree [rad] = pi/180
alias arcmin = 0.00029088820866572158  # arcminute [rad] = degree/60
alias arcminute = 0.00029088820866572158  # arcminute [rad] = degree/60
alias arcsec = 4.8481368110953599e-06  # arcsecond [rad] = arcmin/60
alias arcsecond = 4.8481368110953599e-06  # arcsecond [rad] = arcmin/60

# =======================
# TIME IN SECONDS
# =======================

alias minute = 60.0  # minute [s]
alias hour = 3600.0  # hour [s] = 60*minute
alias day = 86400.0  # day [s] = 24*hour
alias week = 604800.0  # week [s] = 7*day
alias year = 31536000.0  # year [s] = 365*day
alias Julian_year = 31557600.0  # Julian year [s] = 365.25*day

# =======================
# LENGTH IN METERS
# =======================

# Basic length units
alias inch = 0.0254  # inch [m]
alias foot = 0.3048  # foot [m] = 12*inch
alias yard = 0.9144  # yard [m] = 3*foot
alias mile = 1609.344  # mile [m] = 1760*yard
alias mil = 2.54e-05  # mil [m] = inch/1000
alias pt = 0.00035277777777777776  # point [m] = inch/72
alias point = 0.00035277777777777776  # point [m] = inch/72

# Survey units
alias survey_foot = 0.30480060960121924  # survey foot [m] = 1200.0/3937
alias survey_mile = 1609.3472186944375  # survey mile [m] = 5280*survey_foot

# Maritime and scientific units
alias nautical_mile = 1852.0  # nautical mile [m]
alias fermi = 1e-15  # fermi [m]
alias angstrom = 1e-10  # angstrom [m]
alias micron = 1e-06  # micron [m]

# Astronomical units
alias au = 149597870700.0  # astronomical unit [m]
alias astronomical_unit = 149597870700.0  # astronomical unit [m]
alias light_year = 9460730472580800.0  # light year [m] = Julian_year*c
alias parsec = 3.0856775814913673e16  # parsec [m] = au/arcsec

# =======================
# PRESSURE IN PASCALS
# =======================

alias atm = 101325.0  # standard atmosphere [Pa] (2022 CODATA)
alias atmosphere = 101325.0  # standard atmosphere [Pa]
alias bar = 100000.0  # bar [Pa] = 1e5
alias torr = 133.32236842105263  # torr [Pa] = atm/760
alias mmHg = 133.32236842105263  # mmHg [Pa] = atm/760
alias psi = 6894.757293168361  # psi [Pa] = pound*g/(inch*inch)

# =======================
# AREA IN SQUARE METERS
# =======================

alias hectare = 10000.0  # hectare [m^2] = 1e4
alias acre = 4046.8564224  # acre [m^2] = 43560*foot^2

# =======================
# VOLUME IN CUBIC METERS
# =======================

# Metric volume
alias liter = 0.001  # liter [m^3] = 1e-3
alias litre = 0.001  # litre [m^3] = 1e-3

# US volume
alias gallon = 0.003785411784  # gallon (US) [m^3] = 231*inch^3
alias gallon_US = 0.003785411784  # gallon (US) [m^3] = 231*inch^3
alias fluid_ounce = 2.95735295625e-05  # fluid ounce (US) [m^3] = gallon_US/128
alias fluid_ounce_US = 2.95735295625e-05  # fluid ounce (US) [m^3] = gallon_US/128
alias bbl = 0.158987294928  # barrel [m^3] = 42*gallon_US
alias barrel = 0.158987294928  # barrel [m^3] = 42*gallon_US

# UK volume
alias gallon_imp = 0.00454609  # gallon (UK) [m^3] = 4.54609e-3
alias fluid_ounce_imp = 2.8413062499999997e-05  # fluid ounce (UK) [m^3] = gallon_imp/160

# =======================
# SPEED IN METERS PER SECOND
# =======================

alias kmh = 0.2777777777777778  # km/h [m s^-1] = 1e3/hour
alias mph = 0.44704  # mph [m s^-1] = mile/hour
alias mach = 340.5  # Mach [m s^-1] (approx at 15°C, 1 atm)
alias speed_of_sound = 340.5  # speed of sound [m s^-1] (approx at 15°C, 1 atm)
alias knot = 0.5144444444444445  # knot [m s^-1] = nautical_mile/hour

# =======================
# TEMPERATURE IN KELVIN
# =======================

alias zero_Celsius = 273.15  # zero Celsius [K]
alias degree_Fahrenheit = 0.5555555555555556  # degree Fahrenheit [K] (for differences only) = 1/1.8

# =======================
# ENERGY IN JOULES
# =======================

# Basic energy units
alias eV = 1.602176634e-19  # electron volt [J] = elementary_charge
alias electron_volt = 1.602176634e-19  # electron volt [J] = elementary_charge

# Thermal energy units
alias calorie = 4.184  # calorie (thermochemical) [J]
alias calorie_th = 4.184  # calorie (thermochemical) [J]
alias calorie_IT = 4.1868  # calorie (International Steam Table) [J]
alias erg = 1e-07  # erg [J] = 1e-7

# British thermal units
alias Btu_th = 1054.3502644888888  # BTU (thermochemical) [J] = pound*degree_Fahrenheit*calorie_th/gram
alias Btu = 1055.05585262  # BTU (International Steam Table) [J] = pound*degree_Fahrenheit*calorie_IT/gram
alias Btu_IT = 1055.05585262  # BTU (International Steam Table) [J]

# Explosive energy
alias ton_TNT = 4184000000.0  # ton of TNT [J] = 1e9*calorie_th

# =======================
# POWER IN WATTS
# =======================

alias hp = 745.6998715822702  # horsepower [W] = 550*foot*pound*g
alias horsepower = 745.6998715822702  # horsepower [W] = 550*foot*pound*g

# =======================
# FORCE IN NEWTONS
# =======================

alias dyn = 1e-05  # dyne [N] = 1e-5
alias dyne = 1e-05  # dyne [N] = 1e-5
alias lbf = 4.4482216152605  # pound force [N] = pound*g
alias pound_force = 4.4482216152605  # pound force [N] = pound*g
alias kgf = 9.80665  # kilogram force [N] = g
alias kilogram_force = 9.80665  # kilogram force [N] = g
