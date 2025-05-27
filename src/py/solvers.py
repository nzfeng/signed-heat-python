import numpy as np
import shm3d_bindings as shm3db

class SignedHeatTetSolver():

	def __init__(self, verbose):
		self.bound_solver = shm3db.SignedHeatTetSolver(verbose)

	def compute_distance_to_mesh(V, F, level_set_constraint="ZeroSet", t_coef=1., h_coef=0., rebuild=True, scale=2.):
		self.bound_solver.compute_distance_to_mesh(V, F, level_set_constraint, t_coef, h_coef, rebuild, scale)

	def compute_distance_to_point_cloud(P, N, t_coef=1., h_coef=0., rebuild=True, scale=2.):
		self.bound_solver.compute_distance_to_point_cloud(P, N, t_coef, h_coef, rebuild, scale)

class SignedHeatGridSolver():

	def __init__(self, verbose):
		self.bound_solver = shm3db.SignedHeatGridSolver(verbose)

	def compute_distance_to_mesh(V, F, t_coef=1., h_coef=0., rebuild=True, scale=2.):
		self.bound_solver.compute_distance_to_mesh(V, F, level_set_constraint, t_coef, h_coef, rebuild, scale)

	def compute_distance_to_point_cloud(P, N, t_coef=1., h_coef=0., rebuild=True, scale=2.):
		self.bound_solver.compute_distance_to_point_cloud(P, N, t_coef, h_coef, rebuild, scale)