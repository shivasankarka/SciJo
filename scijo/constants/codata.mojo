"""
CODATA Physical Constants Module

This module provides access to the CODATA 2022 recommended values of fundamental
physical constants as published by the Committee on Data for Science and Technology.

The constants are organized in a global dictionary structure following the same
format as SciPy's constants module, ensuring compatibility and familiarity for
scientific computing applications.

Author: Shivasankar K.A
Version: 0.1.0
Date: July 2025

Usage:
    The constants can be accessed through the global dictionary and used in
    scientific calculations requiring precise physical constants.

References:
    - CODATA 2022 Internationally Recommended Values
    - https://github.com/scipy/scipy/blob/main/scipy/constants/_codata.py

Note:
    This implementation is based on the official CODATA 2022 adjustment of
    fundamental physical constants, ensuring the highest accuracy for
    scientific computations.
"""
# =======================
# CODATA PHYSICAL CONSTANTS
# =======================
# Global dictionary of CODATA 2022 physical constants
# Based on https://github.com/scipy/scipy/blob/main/scipy/constants/_codata.py


struct PhysicalConstant[dtype: DType](Copyable, Movable):
    """Physical constant data containing value, unit, and uncertainty."""

    var value: Scalar[dtype]
    var unit: String
    var uncertainty: Scalar[dtype]

    constrained[
        dtype.is_floating_point(),
        "The physical constant must be a floating-point type.",
    ]()

    fn __init__(
        out self, value: Scalar[dtype], unit: String, uncertainty: Scalar[dtype]
    ):
        self.value = value
        self.unit = unit
        self.uncertainty = uncertainty

    fn __str__(self) raises -> String:
        return String("{} {} ± {}").format(
            self.value, self.unit, self.uncertainty
        )


@parameter
fn create_physical_constants_dict[
    dtype: DType
]() -> Dict[String, PhysicalConstant[dtype]]:
    """Create the comprehensive CODATA 2022 physical constants dictionary."""
    var constants = Dict[String, PhysicalConstant[dtype]]()

    # Fundamental physical constants (2022 CODATA values)
    constants["speed_of_light_in_vacuum"] = PhysicalConstant[dtype](
        299792458.0, "m s^-1", 0.0
    )
    constants["elementary_charge"] = PhysicalConstant[dtype](
        1.602176634e-19, "C", 0.0
    )
    constants["Planck_constant"] = PhysicalConstant[dtype](
        6.62607015e-34, "J Hz^-1", 0.0
    )
    constants["reduced_Planck_constant"] = PhysicalConstant[dtype](
        1.0545718176461565e-34, "J s", 0.0
    )
    constants["Boltzmann_constant"] = PhysicalConstant[dtype](
        1.380649e-23, "J K^-1", 0.0
    )
    constants["Avogadro_constant"] = PhysicalConstant[dtype](
        6.02214076e23, "mol^-1", 0.0
    )
    constants["molar_gas_constant"] = PhysicalConstant[dtype](
        8.314462618, "J mol^-1 K^-1", 0.0
    )

    # Electromagnetic constants
    constants["vacuum_electric_permittivity"] = PhysicalConstant[dtype](
        8.8541878188e-12, "F m^-1", 1.3e-21
    )
    constants["vacuum_mag_permeability"] = PhysicalConstant[dtype](
        1.25663706127e-06, "N A^-2", 1.9e-16
    )
    constants["characteristic_impedance_of_vacuum"] = PhysicalConstant[dtype](
        376.730313412, "ohm", 5.9e-8
    )

    # Gravitational constants
    constants["Newtonian_constant_of_gravitation"] = PhysicalConstant[dtype](
        6.6743e-11, "m^3 kg^-1 s^-2", 1.5e-15
    )
    constants["standard_acceleration_of_gravity"] = PhysicalConstant[dtype](
        9.80665, "m s^-2", 0.0
    )

    # Fine structure constants
    constants["fine_structure_constant"] = PhysicalConstant[dtype](
        0.0072973525643, "", 1.1e-12
    )
    constants["inverse_fine_structure_constant"] = PhysicalConstant[dtype](
        137.035999177, "", 2.1e-8
    )

    # Particle masses and properties
    constants["electron_mass"] = PhysicalConstant[dtype](
        9.1093837139e-31, "kg", 2.8e-40
    )
    constants["proton_mass"] = PhysicalConstant[dtype](
        1.67262192595e-27, "kg", 5.2e-37
    )
    constants["neutron_mass"] = PhysicalConstant[dtype](
        1.67492750056e-27, "kg", 8.5e-37
    )
    constants["atomic_mass_constant"] = PhysicalConstant[dtype](
        1.66053906892e-27, "kg", 5.2e-37
    )
    constants["unified_atomic_mass_unit"] = PhysicalConstant[dtype](
        1.66053906892e-27, "kg", 5.2e-37
    )

    # Additional particle masses (from CODATA 2022)
    constants["muon_mass"] = PhysicalConstant[dtype](
        1.883531627e-28, "kg", 4.2e-37
    )
    constants["tau_mass"] = PhysicalConstant[dtype](3.16754e-27, "kg", 2.1e-31)
    constants["alpha_particle_mass"] = PhysicalConstant[dtype](
        6.6446573450e-27, "kg", 2.1e-36
    )
    constants["deuteron_mass"] = PhysicalConstant[dtype](
        3.3435837768e-27, "kg", 1.0e-36
    )
    constants["helion_mass"] = PhysicalConstant[dtype](
        5.0064127862e-27, "kg", 1.6e-36
    )
    constants["triton_mass"] = PhysicalConstant[dtype](
        5.0073567512e-27, "kg", 1.6e-36
    )

    # Atomic units and radii
    constants["Bohr_radius"] = PhysicalConstant[dtype](
        5.29177210544e-11, "m", 8.2e-21
    )
    constants["classical_electron_radius"] = PhysicalConstant[dtype](
        2.8179403205e-15, "m", 1.3e-24
    )
    constants["Thomson_cross_section"] = PhysicalConstant[dtype](
        6.6524587051e-29, "m^2", 6.2e-38
    )

    # Additional atomic units (from CODATA 2022)
    constants["atomic_unit_of_action"] = PhysicalConstant[dtype](
        1.0545718176461565e-34, "J s", 0.0  # exact
    )
    constants["atomic_unit_of_charge"] = PhysicalConstant[dtype](
        1.602176634e-19, "C", 0.0  # exact
    )
    constants["atomic_unit_of_mass"] = PhysicalConstant[dtype](
        9.1093837139e-31, "kg", 2.8e-40
    )
    constants["atomic_unit_of_current"] = PhysicalConstant[dtype](
        6.6236182375082e-3, "A", 7.2e-18
    )
    constants["atomic_unit_of_electric_potential"] = PhysicalConstant[dtype](
        27.211386245981, "V", 3.0e-11
    )
    constants["atomic_unit_of_electric_field"] = PhysicalConstant[dtype](
        5.14220675112e11, "V m^-1", 8.0e4
    )
    constants["atomic_unit_of_electric_dipole_moment"] = PhysicalConstant[
        dtype
    ](8.4783536198e-30, "C m", 1.3e-39)
    constants["atomic_unit_of_electric_polarizability"] = PhysicalConstant[
        dtype
    ](1.64877727212e-41, "C^2 m^2 J^-1", 5.1e-52)
    constants["atomic_unit_of_force"] = PhysicalConstant[dtype](
        8.2387235038e-8, "N", 1.3e-15
    )
    constants["atomic_unit_of_velocity"] = PhysicalConstant[dtype](
        2.18769126216e6, "m s^-1", 3.4e-1
    )
    constants["atomic_unit_of_momentum"] = PhysicalConstant[dtype](
        1.99285191545e-24, "kg m s^-1", 3.1e-34
    )

    # Magnetic moments and related constants
    constants["Bohr_magneton"] = PhysicalConstant[dtype](
        9.2740100657e-24, "J T^-1", 2.9e-33
    )
    constants["nuclear_magneton"] = PhysicalConstant[dtype](
        5.0507837393e-27, "J T^-1", 1.6e-36
    )
    constants["electron_mag_mom"] = PhysicalConstant[dtype](
        -9.2847646917e-24, "J T^-1", 2.9e-33
    )
    constants["proton_mag_mom"] = PhysicalConstant[dtype](
        1.41060679545e-26, "J T^-1", 6.0e-36
    )

    # Additional magnetic moments (from CODATA 2022)
    constants["neutron_magnetic_moment"] = PhysicalConstant[dtype](
        -9.6623653e-27, "J T^-1", 2.3e-33
    )
    constants["muon_magnetic_moment"] = PhysicalConstant[dtype](
        -4.49044830e-26, "J T^-1", 1.0e-33
    )
    constants["deuteron_magnetic_moment"] = PhysicalConstant[dtype](
        4.330735087e-27, "J T^-1", 1.1e-35
    )
    constants["helion_magnetic_moment"] = PhysicalConstant[dtype](
        -1.07461755198e-26, "J T^-1", 9.3e-35
    )
    constants["triton_magnetic_moment"] = PhysicalConstant[dtype](
        1.5046095178e-26, "J T^-1", 3.0e-35
    )

    # Gyromagnetic ratios
    constants["electron_gyromagnetic_ratio"] = PhysicalConstant[dtype](
        1.76085962784e11, "s^-1 T^-1", 5.5e4
    )
    constants["proton_gyromagnetic_ratio"] = PhysicalConstant[dtype](
        2.6752218708e8, "s^-1 T^-1", 1.1e1
    )
    constants["neutron_gyromagnetic_ratio"] = PhysicalConstant[dtype](
        1.83247174e8, "s^-1 T^-1", 4.3e1
    )

    # Atomic and molecular constants
    constants["Rydberg_constant"] = PhysicalConstant[dtype](
        10973731.568157, "m^-1", 1.2e-5
    )
    constants["Hartree_energy"] = PhysicalConstant[dtype](
        4.3597447222060e-18, "J", 4.8e-28
    )
    constants["atomic_unit_of_length"] = PhysicalConstant[dtype](
        5.29177210544e-11, "m", 8.2e-21
    )
    constants["atomic_unit_of_energy"] = PhysicalConstant[dtype](
        4.3597447222060e-18, "J", 4.8e-28
    )
    constants["atomic_unit_of_time"] = PhysicalConstant[dtype](
        2.4188843265864e-17, "s", 2.6e-27
    )

    # Compton wavelengths
    constants["compton_wavelength"] = PhysicalConstant[dtype](
        2.42631023538e-12, "m", 7.6e-22
    )
    constants["electron_compton_wavelength"] = PhysicalConstant[dtype](
        2.42631023538e-12, "m", 7.6e-22
    )
    constants["proton_compton_wavelength"] = PhysicalConstant[dtype](
        1.32140985360e-15, "m", 4.1e-25
    )
    constants["neutron_compton_wavelength"] = PhysicalConstant[dtype](
        1.31959090382e-15, "m", 6.7e-25
    )

    # Thermodynamic constants
    constants["stefan_boltzmann_constant"] = PhysicalConstant[dtype](
        5.670374419e-08, "W m^-2 K^-4", 0.0
    )
    constants["wien_wavelength_displacement_law_constant"] = PhysicalConstant[
        dtype
    ](0.002897771955, "m K", 0.0)
    constants["wien_frequency_displacement_law_constant"] = PhysicalConstant[
        dtype
    ](5.878925757e10, "Hz K^-1", 0.0)

    # Radiation constants
    constants["first_radiation_constant"] = PhysicalConstant[dtype](
        3.741771852e-16, "W m^2", 0.0
    )
    constants["second_radiation_constant"] = PhysicalConstant[dtype](
        0.01438776877, "m K", 0.0
    )

    # Quantum electrical constants
    constants["conductance_quantum"] = PhysicalConstant[dtype](
        7.748091729e-05, "S", 0.0
    )
    constants["inverse_of_conductance_quantum"] = PhysicalConstant[dtype](
        12906.40372, "ohm", 0.0
    )
    constants["magnetic_flux_quantum"] = PhysicalConstant[dtype](
        2.067833848e-15, "Wb", 0.0
    )
    constants["Josephson_constant"] = PhysicalConstant[dtype](
        4.835978484e14, "Hz V^-1", 0.0
    )
    constants["von_klitzing_constant"] = PhysicalConstant[dtype](
        25812.80745, "ohm", 0.0
    )

    # Additional quantum constants (from CODATA 2022)
    constants["quantum_of_circulation"] = PhysicalConstant[dtype](
        3.6369475467e-4, "m^2 s^-1", 1.1e-13
    )
    constants["Faraday_constant"] = PhysicalConstant[dtype](
        96485.33212, "C mol^-1", 0.0  # exact
    )

    # Natural units
    constants["natural_unit_of_length"] = PhysicalConstant[dtype](
        3.8615926744e-13, "m", 1.2e-22
    )
    constants["natural_unit_of_energy"] = PhysicalConstant[dtype](
        8.1871057880e-14, "J", 2.6e-23
    )
    constants["natural_unit_of_time"] = PhysicalConstant[dtype](
        1.28808866644e-21, "s", 4.0e-31
    )
    constants["natural_unit_of_velocity"] = PhysicalConstant[dtype](
        299792458.0, "m s^-1", 0.0  # exact (speed of light)
    )

    # Mass ratios and conversion factors
    constants["proton_electron_mass_ratio"] = PhysicalConstant[dtype](
        1836.152673426, "", 3.2e-8
    )
    constants["neutron_electron_mass_ratio"] = PhysicalConstant[dtype](
        1838.68366200, "", 7.4e-7
    )
    constants["deuteron_electron_mass_ratio"] = PhysicalConstant[dtype](
        3670.482967655, "", 6.3e-8
    )
    constants["alpha_particle_electron_mass_ratio"] = PhysicalConstant[dtype](
        7294.29954144, "", 2.4e-7
    )

    # Additional mass ratios (from CODATA 2022)
    constants["muon_electron_mass_ratio"] = PhysicalConstant[dtype](
        206.7682827, "", 4.6e-6
    )
    constants["tau_electron_mass_ratio"] = PhysicalConstant[dtype](
        3477.23, "", 2.3e-1
    )
    constants["proton_neutron_mass_ratio"] = PhysicalConstant[dtype](
        0.99862347797, "", 4.0e-10
    )
    constants["deuteron_proton_mass_ratio"] = PhysicalConstant[dtype](
        1.9990075012699, "", 8.4e-12
    )
    constants["helion_proton_mass_ratio"] = PhysicalConstant[dtype](
        2.993152671552, "", 7.0e-11
    )

    # Energy conversion factors
    constants["electron_volt"] = PhysicalConstant[dtype](
        1.602176634e-19, "J", 0.0
    )
    constants["atomic_mass_unit_electron_volt_relationship"] = PhysicalConstant[
        dtype
    ](9.3149410372e8, "eV", 2.9e-2)
    constants["electron_volt_joule_relationship"] = PhysicalConstant[dtype](
        1.602176634e-19, "J", 0.0
    )
    constants["electron_volt_kelvin_relationship"] = PhysicalConstant[dtype](
        1.160451812e4, "K", 0.0
    )

    # Additional conversion factors (from CODATA 2022)
    constants["hartree_electron_volt_relationship"] = PhysicalConstant[dtype](
        27.211386245981, "eV", 3.0e-11
    )
    constants["hartree_joule_relationship"] = PhysicalConstant[dtype](
        4.3597447222060e-18, "J", 4.8e-28
    )
    constants["electron_volt_hertz_relationship"] = PhysicalConstant[dtype](
        2.417989242e14, "Hz", 0.0  # exact
    )
    constants["electron_volt_inverse_meter_relationship"] = PhysicalConstant[
        dtype
    ](
        8.065543937e5, "m^-1", 0.0  # exact
    )

    # Planck units
    constants["planck_length"] = PhysicalConstant[dtype](
        1.616255e-35, "m", 1.8e-40
    )
    constants["planck_mass"] = PhysicalConstant[dtype](
        2.176434e-8, "kg", 2.4e-13
    )
    constants["planck_time"] = PhysicalConstant[dtype](
        5.391247e-44, "s", 6.0e-49
    )
    constants["planck_temperature"] = PhysicalConstant[dtype](
        1.416784e32, "K", 1.6e27
    )

    # Pressure and atmosphere
    constants["standard_atmosphere"] = PhysicalConstant[dtype](
        101325.0, "Pa", 0.0
    )
    constants["standard_state_pressure"] = PhysicalConstant[dtype](
        100000.0, "Pa", 0.0
    )

    # Molar properties and additional constants
    constants["molar_mass_constant"] = PhysicalConstant[dtype](
        1.00000000105e-3, "kg mol^-1", 3.1e-13
    )
    constants["molar_mass_of_carbon_12"] = PhysicalConstant[dtype](
        1.20000000126e-2, "kg mol^-1", 3.7e-12
    )
    constants[
        "molar_volume_of_ideal_gas_(273.15 K, 100 kPa)"
    ] = PhysicalConstant[dtype](2.271095464e-2, "m^3 mol^-1", 0.0)
    constants[
        "molar_volume_of_ideal_gas_(273.15 K, 101.325 kPa)"
    ] = PhysicalConstant[dtype](2.241396954e-2, "m^3 mol^-1", 0.0)

    # Additional fundamental constants (from CODATA 2022)
    constants["hyperfine_transition_frequency_of_Cs_133"] = PhysicalConstant[
        dtype
    ](
        9192631770.0, "Hz", 0.0  # exact
    )
    constants["luminous_efficacy"] = PhysicalConstant[dtype](
        683.0, "lm W^-1", 0.0  # exact
    )
    constants["Loschmidt_constant_273_15_K_100_kPa"] = PhysicalConstant[dtype](
        2.651645804e25, "m^-3", 0.0  # exact
    )
    constants["Loschmidt_constant_273_15_K_101_325_kPa"] = PhysicalConstant[
        dtype
    ](
        2.686780111e25, "m^-3", 0.0  # exact
    )

    # Additional Boltzmann constant relationships
    constants["Boltzmann_constant_in_eV_K"] = PhysicalConstant[dtype](
        8.617333262e-5, "eV K^-1", 0.0  # exact
    )
    constants["Boltzmann_constant_in_Hz_K"] = PhysicalConstant[dtype](
        2.083661912e10, "Hz K^-1", 0.0  # exact
    )

    # Lattice and crystallographic constants
    constants["lattice_parameter_of_silicon"] = PhysicalConstant[dtype](
        5.431020511e-10, "m", 8.9e-18
    )
    constants["lattice_spacing_of_ideal_Si_220"] = PhysicalConstant[dtype](
        1.920155716e-10, "m", 3.2e-18
    )

    return constants


# Compile time dictionary of physical constants
alias physical_constants = create_physical_constants_dict[DType.float64]()


fn value(key: String) raises -> Scalar[DType.float64]:
    """
    Get the value of a physical constant.

    Args:
        key: Name of the physical constant.

    Returns:
        The numerical value of the constant.
    """
    if key in physical_constants:
        return physical_constants[key].value
    else:
        print("Warning: Unknown constant '" + key + "'")
        return 0.0


fn unit(key: String) raises -> String:
    """
    Get the unit of a physical constant.

    Args:
        key: Name of the physical constant.

    Returns:
        The unit string of the constant.
    """
    if key in physical_constants:
        return physical_constants[key].unit
    else:
        print("Warning: Unknown constant '" + key + "'")
        return ""


fn precision(key: String) raises -> Scalar[DType.float64]:
    """
    Get the relative precision (uncertainty/value) of a physical constant.

    Args:
        key: Name of the physical constant.

    Returns:
        The relative precision of the constant.
    """
    if key in physical_constants:
        var constant = physical_constants[key]
        if constant.value != 0.0:
            return constant.uncertainty / constant.value
        else:
            return 0.0
    else:
        print("Warning: Unknown constant '" + key + "'")
        return 0.0


fn find(substring: String = "") raises -> List[String]:
    """
    Find physical constants containing a substring in their name.

    Args:
        substring: Substring to search for (empty returns all constants).

    Returns:
        List of constant names containing the substring.
    """
    var result = List[String]()

    for item in physical_constants.items():
        var key = item.key
        if substring == "" or substring in key:
            result.append(key)

    return result


# Additional helper functions for common access patterns
fn get_constant_tuple(
    key: String,
) raises -> Tuple[Scalar[DType.float64], String, Scalar[DType.float64]]:
    """
    Get a physical constant as a tuple (value, unit, uncertainty).

    Args:
        key: Name of the physical constant.

    Returns:
        Tuple containing (value, unit, uncertainty).
    """
    if key in physical_constants:
        var constant = physical_constants[key]
        return (constant.value, constant.unit, constant.uncertainty)
    else:
        print("Warning: Unknown constant '" + key + "'")
        return (0.0, "", 0.0)


fn list_all_constants() raises -> List[String]:
    """
    Get a list of all available physical constant names.

    Returns:
        List of all constant names in the database.
    """
    var result = List[String]()
    for item in physical_constants.items():
        result.append(item.key)
    return result
