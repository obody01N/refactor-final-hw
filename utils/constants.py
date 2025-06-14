amu = 1.66053906660e-27  # kg, atomic mass unit
angstrom = 1e-10
eV = 1.602176634e-19  # Joules
k_B = 1.380649e-23  # J/K
mass=12.0
mass_kg = mass * amu  # kg

dt = 1e-14  # seconds
v_unit = 1/(angstrom / dt) # = 1e-4
a_unit = eV/(angstrom*amu) / (angstrom / dt**2) # = 0.9648533215665326