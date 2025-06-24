import yaml
from pathlib import Path
import numpy as np

class Config:
    def __init__(self, config_file="config_file.yaml"):
        with open(config_file, 'r') as f:
            self.data = yaml.safe_load(f)
        
        # Simulation parameters
        sim = self.data['simulation']
        self.nx, self.ny, self.nz = map(int, sim['dimensions'])
        self.dt = float(sim['timestep'])
        self.steps = int(sim['steps'])
        self.temperature = list(map(float, sim['temperature']))
        self.mass = float(sim['mass'])
        self.boundary_condition = sim['boundary_condition']
        self.track_particle = bool(sim['track_particle'])
        self.write_output = bool(sim['write_output'])
        
        # Potential parameters
        pot = self.data['potentials']
        self.two_potentials = bool(pot['two_potentials'])
        self.epsilon = float(pot['lj']['epsilon'])
        self.sigma = float(pot['lj']['sigma'])
        self.cutoff = float(pot['lj']['cutoff'])
        self.ro = float(pot['lj']['smooth_ro'])
        self.rc = float(pot['lj']['smooth_rc'])
        self.D_e = float(pot['morse']['D_e'])
        self.a = float(pot['morse']['a'])
        self.r_e = float(pot['morse']['r_e'])
        self.cutoff_face = float(pot['cutoff_face'])
        
        # Lattice parameters
        lat = self.data['lattice']
        self.basis_vectors = np.array(lat['basis_vectors'], dtype=float)
        self.basis_atoms = np.array(lat['basis_atoms'], dtype=float)
        
        # Analysis parameters
        ana = self.data['analysis']
        self.r_max = float(ana['rdf']['r_max'])
        self.dr = float(ana['rdf']['dr'])
        self.T_start = float(ana['heat_capacity']['T_start'])
        self.T_end = float(ana['heat_capacity']['T_end'])
        self.T_number = int(ana['heat_capacity']['T_number'])

    def save(self, filename="config_saved.yaml"):
        with open(filename, 'w') as f:
            yaml.dump(self.data, f)
