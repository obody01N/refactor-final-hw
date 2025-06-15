import yaml
from pathlib import Path
import numpy as np

class Config:
    def __init__(self, config_file="config.yaml"):
        with open(config_file, 'r') as f:
            self.data = yaml.safe_load(f)
        
        # Simulation parameters
        sim = self.data['simulation']
        self.nx, self.ny, self.nz = sim['dimensions']
        self.dt = sim['timestep']
        self.steps = sim['steps']
        self.temperature = sim['temperature']
        self.mass = sim['mass']
        self.boundary_condition = sim['boundary_condition']
        self.track_particle = sim['track_particle']
        self.write_output = sim['write_output']
        
        # Potential parameters
        pot = self.data['potentials']
        self.epsilon = pot['lj']['epsilon']
        self.sigma = pot['lj']['sigma']
        self.cutoff = pot['lj']['cutoff']
        self.ro = pot['lj']['smooth_ro']
        self.rc = pot['lj']['smooth_rc']
        self.D_e = pot['morse']['D_e']
        self.a = pot['morse']['a']
        self.r_e = pot['morse']['r_e']
        self.cutoff_face = pot['cutoff_face']
        
        # Lattice parameters
        lat = self.data['lattice']
        self.basis_vectors = np.array(lat['basis_vectors'])
        self.basis_atoms = np.array(lat['basis_atoms'])
        
        # Analysis parameters
        ana = self.data['analysis']
        self.r_max = ana['rdf']['r_max']
        self.dr = ana['rdf']['dr']
        self.T_start = ana['heat_capacity']['T_start']
        self.T_end = ana['heat_capacity']['T_end']
        self.T_number = ana['heat_capacity']['T_number']

    def save(self, filename="config_saved.yaml"):
        with open(filename, 'w') as f:
            yaml.dump(self.data, f)

if __name__ == "__main__":
    config = Config()
    print("Configuration loaded successfully.")
    print(f"Simulation steps: {config.steps}")
    print(f"Temperature range: {config.T_start} K to {config.T_end} K")
    print(f"Lattice basis vectors: {config.basis_vectors}")
    print(f"Potential epsilon: {config.epsilon}, sigma: {config.sigma}")
    config.save("config_saved.yaml")  # Save the configuration to a file if needed