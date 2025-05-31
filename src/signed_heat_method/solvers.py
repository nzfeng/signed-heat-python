import numpy as np
import shm3d_bindings as shm3db

class SignedHeatTetSolver():

	def __init__(self, verbose=True) -> None:
		self.bound_solver = shm3db.SignedHeatTetSolver(verbose)

	def get_vertices(self) -> np.ndarray:
		return self.bound_solver.get_vertices()

	def get_tets(self) -> np.ndarray:
		return self.bound_solver.get_tets()

	def compute_distance_to_mesh(self, V: np.ndarray, F: list[list[int]], level_set_constraint: str="ZeroSet", t_coef: float=1., h_coef: float=0., rebuild: bool=True) -> np.ndarray:
		return self.bound_solver.compute_distance_to_mesh(V, F, level_set_constraint, t_coef, h_coef, rebuild)

	def compute_distance_to_point_cloud(self, P: np.ndarray, N: np.ndarray, level_set_constraint: str="ZeroSet", t_coef: float=1., h_coef: float=0., rebuild: bool=True) -> np.ndarray:
		return self.bound_solver.compute_distance_to_point_cloud(P, N, level_set_constraint, t_coef, h_coef, rebuild)

	def isosurface(self, phi: np.ndarray, isoval: float=0.): # TODO: -> tuple(np.ndarray, list[list[int]])
		return self.bound_solver.isosurface(phi, isoval)

class SignedHeatGridSolver():

	def __init__(self, verbose: bool=True) -> None:
		self.bound_solver = shm3db.SignedHeatGridSolver(verbose)

	def get_grid_resolution(self) -> list[int]:
		return self.bound_solver.get_grid_resolution()

	def get_bbox(self): # TODO: -> tuple(np.ndarray, np.ndarray)
		return self.bound_solver.get_bbox()

	def to_grid_array(self, phi: np.ndarray) -> np.ndarray:
		'''
		Convert an array of size (dim_x * dim_y * dim_z) to a NumPy array of shape (dim_x, dim_y, dim_z).
		Warning: Logic is duplicated between here and "indicesToNodeIndex()" in signed-heat-3d/src/signed_heat_grid_solver.cpp.
		'''
		nx, ny, nz = self.get_grid_resolution()
		i_idx, j_idx, k_idx = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')
		indices = i_idx + j_idx * ny + k_idx * (nx * ny)
		new_array = phi[indices]
		return new_array

	def compute_distance_to_mesh(self, V: np.ndarray, F: list[list[int]], t_coef: float=1., h_coef: float=0., rebuild: bool=True) -> np.ndarray:
		return self.bound_solver.compute_distance_to_mesh(V, F, t_coef, h_coef, rebuild)

	def compute_distance_to_point_cloud(self, P: np.ndarray, N: np.ndarray, t_coef: float=1., h_coef: float=0., rebuild: bool=True) -> np.ndarray:
		return self.bound_solver.compute_distance_to_point_cloud(P, N, t_coef, h_coef, rebuild)