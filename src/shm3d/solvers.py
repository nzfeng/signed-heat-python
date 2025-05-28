import numpy as np
import shm3d_bindings as shm3db

class SignedHeatTetSolver():

	def __init__(self, verbose):
		self.bound_solver = shm3db.SignedHeatTetSolver(verbose)

	def get_vertices(self):
		return self.bound_solver.get_vertices()

	def get_tets(self):
		return self.bound_solver.get_tets()

	def compute_distance_to_mesh(self, V, F, level_set_constraint="ZeroSet", t_coef=1., h_coef=0., rebuild=True):
		return self.bound_solver.compute_distance_to_mesh(V, F, level_set_constraint, t_coef, h_coef, rebuild)

	def compute_distance_to_point_cloud(self, P, N, level_set_constraint="ZeroSet", t_coef=1., h_coef=0., rebuild=True):
		return self.bound_solver.compute_distance_to_point_cloud(P, N, level_set_constraint, t_coef, h_coef, rebuild)

	def isosurface(self, phi, isoval):
		return self.bound_solver.isosurface(phi, isoval)

class SignedHeatGridSolver():

	def __init__(self, verbose):
		self.bound_solver = shm3db.SignedHeatGridSolver(verbose)

	def get_grid_resolution(self):
		return self.bound_solver.get_grid_resolution()

	def get_bbox(self):
		return self.bound_solver.get_bbox()

	def to_grid_array(self, phi):
		'''
		Convert an array of size (dim_x * dim_y * dim_z) to a NumPy array of shape (dim_x, dim_y, dim_z).
		Warning: Logic is duplicated between here and "indicesToNodeIndex()" in signed-heat-3d/src/signed_heat_grid_solver.cpp!
		'''
		nx, ny, nz = self.get_grid_resolution()
		new_array = np.empty((nx, ny, nz))
		for i in range(nx):
			for j in range(ny):
				for k in range(nz):
					idx = i + j * ny + k * (nx * ny)
					new_array[i,j,k] = phi[idx]

		return new_array

	def compute_distance_to_mesh(self, V, F, t_coef=1., h_coef=0., rebuild=True):
		return self.bound_solver.compute_distance_to_mesh(V, F, t_coef, h_coef, rebuild)

	def compute_distance_to_point_cloud(self, P, N, t_coef=1., h_coef=0., rebuild=True):
		return self.bound_solver.compute_distance_to_point_cloud(P, N, t_coef, h_coef, rebuild)