import numpy as np
import shm3d_bindings as shm3db

# Warning: Default values are duplicated between here and `SignedHeat3DOptions` in signed-heat-3d/include/signed_heat_3d.h


class SignedHeatTetSolver:
	def __init__(self, verbose=True) -> None:
		self.bound_solver = shm3db.SignedHeatTetSolver(verbose)

	def get_vertices(self) -> np.ndarray:
		return self.bound_solver.get_vertices()

	def get_tets(self) -> np.ndarray:
		return self.bound_solver.get_tets()

	def compute_distance_to_mesh(
		self,
		V: np.ndarray,
		F: list[list[int]],
		options: dict = {},
	) -> np.ndarray:
		return self.bound_solver.compute_distance_to_mesh(
			V,
			F,
			options.get('level_set_constraint', 'ZeroSet'),
			options.get('t_coef', 1.0),
			options.get('bbox_min', np.array([1.0, 1.0, 1.0])),
			options.get('bbox_max', np.array([-1.0, -1.0, -1.0])),
			options.get('resolution', np.array([0, 0, 0])),
		)

	def compute_distance_to_point_cloud(self, P: np.ndarray, N: np.ndarray, options: dict = {}) -> np.ndarray:
		return self.bound_solver.compute_distance_to_point_cloud(
			P,
			N,
			options.get('level_set_constraint', 'ZeroSet'),
			options.get('t_coef', 1.0),
			options.get('bbox_min', np.array([1.0, 1.0, 1.0])),
			options.get('bbox_max', np.array([-1.0, -1.0, -1.0])),
			options.get('resolution', np.array([0, 0, 0])),
		)

	def isosurface(self, phi: np.ndarray, isoval: float = 0.0) -> tuple[np.ndarray, list[list[int]]]:
		return self.bound_solver.isosurface(phi, isoval)


class SignedHeatGridSolver:
	def __init__(self, verbose: bool = True) -> None:
		self.bound_solver = shm3db.SignedHeatGridSolver(verbose)

	def get_grid_resolution(self) -> list[int]:
		return self.bound_solver.get_grid_resolution()

	def get_bbox(self) -> tuple[np.ndarray, np.ndarray]:
		return self.bound_solver.get_bbox()

	def to_grid_array(self, phi: np.ndarray) -> np.ndarray:
		"""
		Convert an array of size (dim_x * dim_y * dim_z) to a NumPy array of shape (dim_x, dim_y, dim_z).
		Warning: Logic is duplicated between here and "indicesToNodeIndex()" in signed-heat-3d/src/signed_heat_grid_solver.cpp.
		"""
		nx, ny, nz = self.get_grid_resolution()
		i_idx, j_idx, k_idx = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')
		indices = i_idx + j_idx * nx + k_idx * (nx * ny)
		new_array = phi[indices]
		return new_array

	def compute_distance_to_mesh(self, V: np.ndarray, F: list[list[int]], options: dict = {}) -> np.ndarray:
		return self.bound_solver.compute_distance_to_mesh(
			V,
			F,
			options.get('t_coef', 1.0),
			options.get('bbox_min', np.array([1.0, 1.0, 1.0])),
			options.get('bbox_max', np.array([-1.0, -1.0, -1.0])),
			options.get('resolution', np.array([0, 0, 0])),
		)

	def compute_distance_to_point_cloud(
		self, P: np.ndarray, N: np.ndarray, t_coef: float = 1.0, options: dict = {}
	) -> np.ndarray:
		return self.bound_solver.compute_distance_to_point_cloud(
			P,
			N,
			options.get('t_coef', 1.0),
			options.get('bbox_min', np.array([1.0, 1.0, 1.0])),
			options.get('bbox_max', np.array([-1.0, -1.0, -1.0])),
			options.get('resolution', np.array([0, 0, 0])),
		)
