import os
import sys
import numpy as np
import potpourri3d as pp3d
import platform

from mesh_io import *

# Path to where the bindings live
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
if os.name == 'nt':  # if Windows
	# handle default location where VS puts binary
	sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build', 'Debug')))
else:
	# normal / unix case
	sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'build')))

asset_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '.', 'assets'))

import signedheat3d as shm

TET_RESOLUTION = np.array([8, 8, 8]) if platform.system != 'linux' else np.array([2, 2, 2])


def area_weighted_vertex_normals(V: np.ndarray, F: list[list[int]]) -> np.ndarray:
	"""
	Compute vertex normals.
	"""
	# Compute unit normals for each polygon face.
	triangles = V[F]
	edge1 = triangles[:, 1] - triangles[:, 0]
	edge2 = triangles[:, 2] - triangles[:, 0]
	face_normals = np.cross(edge1, edge2)  # area-weighted
	# Compute vertex normals.
	vertex_normals = np.zeros_like(V)
	face_normals_repeated = np.repeat(face_normals, [len(f) for f in F], axis=0)
	F_flat = [i for f in F for i in f]
	np.add.at(vertex_normals, F_flat, face_normals_repeated)
	vertex_normals = vertex_normals / np.linalg.norm(vertex_normals, axis=1, keepdims=True)
	return vertex_normals


def approximate_signed_distance(Q: np.ndarray, P: np.ndarray, N: np.ndarray) -> np.ndarray:
	"""
	Approximate signed distance from a set of query points to an oriented point set, using pseudonormal distance.

	Args:
		Q: _ x 3 NumPy array (query points)
		P: _ x 3 NumPy array
		N: _ x 3 NumPy array

	Returns:
		|Q| x 3 NumPy array
	"""
	Q_expanded = Q[:, np.newaxis, :]
	P_expanded = P[np.newaxis, :, :]
	N_expanded = N[np.newaxis, :, :]
	diff = Q_expanded - P_expanded
	all_distances = np.linalg.norm(diff, axis=2)  # (|Q|, |P|)
	dot_products = np.sum(N_expanded * diff, axis=2)  # (|Q|, |P|)
	signed_all_distances = np.sign(dot_products) * all_distances  # (|Q|, |P|)
	# Find closest point by unsigned distance
	abs_distances = np.abs(signed_all_distances)
	min_indices = np.argmin(abs_distances, axis=1)  # (|Q|,)
	signed_distances = signed_all_distances[np.arange(Q.shape[0]), min_indices]
	return signed_distances


def grid_node_positions(nx: int, ny: int, nz: int, bbox_min: np.ndarray, bbox_max: np.ndarray) -> np.ndarray:
	"""
	Return grid positions as a 2D N x 3 NumPy array, where N = nx * ny * nz is the total number of nodes in the grid.
	"""
	x_coords = np.linspace(bbox_min[0], bbox_max[0], nx)
	y_coords = np.linspace(bbox_min[1], bbox_max[1], ny)
	z_coords = np.linspace(bbox_min[2], bbox_max[2], nz)
	i_grid, j_grid, k_grid = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')
	i_flat = i_grid.flatten()
	j_flat = j_grid.flatten()
	k_flat = k_grid.flatten()
	linear_indices = i_flat + j_flat * ny + k_flat * (nx * ny)
	x_positions = x_coords[i_flat]
	y_positions = y_coords[j_flat]
	z_positions = z_coords[k_flat]
	Q = np.zeros((nx * ny * nz, 3))
	Q[linear_indices] = np.column_stack([x_positions, y_positions, z_positions])
	return Q


def test_read_polygon_mesh() -> None:
	V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
	assert len(V.shape) == 2, 'Vertex array should be a 2D NumPy array.'
	assert V.shape[1] == 3, 'Vertex array should be a _ x 3 NumPy array.'
	assert all(len(polygon) == 3 for polygon in F), 'bunny_small should be a triangle mesh.'
	max_elem = max(max(F))
	assert max_elem < V.shape[0], (
		'There is a face with an out-of-bounds vertex index. Faces should be zero-based arrays of indices into vertices.'
	)

	V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'spot.obj'))
	assert len(V.shape) == 2, 'Vertex array should be a 2D NumPy array.'
	assert V.shape[1] == 3, 'Vertex array should be a _ x 3 NumPy array.'
	assert all(len(polygon) == 4 for polygon in F), 'spot should be a quad mesh.'
	max_elem = np.amax(F)
	min_elem = np.amin(F)
	assert max_elem < V.shape[0], (
		'There is a face with an out-of-bounds vertex index. Faces should be zero-based arrays of indices into vertices.'
	)
	assert min_elem >= 0, (
		'There is a face with a negative vertex index. Faces should be zero-based arrays of indices into vertices.'
	)


def test_write_surface_mesh() -> None:
	V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
	out_filepath = os.path.join(asset_path, 'test_mesh.obj')
	write_surface_mesh(V, F, out_filepath)
	V_new, F_new = read_polygon_mesh(os.path.join(asset_path, 'test_mesh.obj'))
	assert np.amax(np.abs(V - V_new)) < 1e-6, 'Vertices of written mesh are different than that of read mesh.'
	assert F == F_new, 'Faces of written mesh are different than that of read mesh.'


def test_read_point_cloud() -> None:
	P, N = read_point_cloud(os.path.join(asset_path, 'bunny.pc'))
	assert len(P.shape) == 2, 'Point array should be a 2D NumPy array.'
	assert P.shape[1] == 3, 'Point array should be a _ x 3 NumPy array.'
	assert len(N.shape) == 2, 'Point normal array should be a 2D NumPy array.'
	assert N.shape[1] == 3, 'Point normal array should be a _ x 3 NumPy array.'
	assert P.shape == N.shape, 'Point and normal arrays should have the same shape.'


class TestTetSolver:
	"""
	The resolution of the tet mesh is set super low because the Ubuntu runner on Github seems it can't handle the memory requirements of TetGen.
	"""

	def test_get_vertices(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatTetSolver(verbose=True)
		solve_options = {'resolution': TET_RESOLUTION}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		vertices = solver.get_vertices()
		assert len(vertices.shape) == 2, 'Vertex array of tet mesh should be a 2D NumPy array.'
		assert vertices.shape[1] == 3, 'Vertex array of tet mesh should be a _ x 3 NumPy array.'
		assert vertices.shape[0] == len(phi), (
			'SDF should have same length as number of vertices in the tet mesh domain.'
		)

	def test_get_tets(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatTetSolver(verbose=True)
		solve_options = {'resolution': TET_RESOLUTION}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		vertices = solver.get_vertices()
		tets = solver.get_tets()
		assert len(tets.shape) == 2, 'Tet array of tet mesh should be a 2D NumPy array.'
		assert tets.shape[1] == 4, 'Tet array of tet mesh should be a _ x 4 NumPy array.'
		max_elem = np.amax(tets)
		min_elem = np.amin(tets)
		assert max_elem < vertices.shape[0], (
			'There is a tet with an out-of-bounds vertex index. Tets should be zero-based arrays of indices into vertices.'
		)
		assert min_elem >= 0, (
			'There is a tet with a negative vertex index. Tets should be zero-based arrays of indices into vertices.'
		)

	def test_compute_distance_to_mesh(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatTetSolver(verbose=True)
		solve_options = {'resolution': TET_RESOLUTION}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		assert len(phi.shape) == 1, 'SDF should be a 1D NumPy array.'
		# Approximate the ground-truth distance, assuming test mesh is closed & perfect.
		N = area_weighted_vertex_normals(V, F)
		# Approximate distance by just using signed distance to vertices (instead of triangle faces.)
		signed_distances = approximate_signed_distance(solver.get_vertices(), V, N)
		span = np.amax(signed_distances) - np.amin(signed_distances)
		residual = (phi - signed_distances) / span
		residual[np.isnan(residual)] = 0.0
		assert np.mean(residual) < 2e-2, 'SDF not close to approximate ground-truth.'

	def test_compute_distance_to_point_cloud(self) -> None:
		P, N = read_point_cloud(os.path.join(asset_path, 'bunny.pc'))
		solver = shm.SignedHeatTetSolver(verbose=True)
		solve_options = {'resolution': TET_RESOLUTION}
		phi = solver.compute_distance_to_point_cloud(P=P, N=N, options=solve_options)
		assert len(phi.shape) == 1, 'SDF should be a 1D NumPy array.'
		# Make sure distance is close to naive distance
		# Approximate signed distance to point cloud using pseudonormal distance.
		signed_distances = approximate_signed_distance(solver.get_vertices(), P, N)
		signed_distances[np.isnan(signed_distances)] = 0.0
		span = np.amax(signed_distances) - np.amin(signed_distances)
		residual = (phi - signed_distances) / span
		residual[np.isnan(residual)] = 0.0
		assert np.mean(residual) < 2e-2, 'SDF not close to approximate ground-truth.'

	def average_squared_distance(
		self, V1: np.ndarray, F1: list[list[int]], V2: np.ndarray, F2: list[list[int]]
	) -> float:
		"""
		Approximate since we'll be using distances to point sets.
		"""
		N1 = area_weighted_vertex_normals(V1, F1)
		N2 = area_weighted_vertex_normals(V2, F2)
		D2_1 = np.square(approximate_signed_distance(V1, V2, N2))
		D2_2 = np.square(approximate_signed_distance(V2, V1, N1))
		return (np.sum(D2_1) + np.sum(D2_2)) / (V1.shape[0] + V2.shape[0])

	def test_isosurface(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatTetSolver(verbose=True)
		solve_options = {'resolution': TET_RESOLUTION, 'level_set_constraint': 'ZeroSet'}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		V_iso, F_iso = solver.isosurface(phi, 0.0)
		error = self.average_squared_distance(V, F, V_iso, F_iso)
		assert error < 1e-2, 'Zero level set is far away from original closed mesh.'


class TestGridSolver:
	def test_get_grid_resolution(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatGridSolver(verbose=True)
		solve_options = {}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		grid_res = solver.get_grid_resolution()
		print(grid_res)
		assert len(grid_res) == 3, 'Grid should be 3D.'
		assert all(res > 0 for res in grid_res), 'Grid should have nonzero length in each dimension.'

	def test_get_bbox(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatGridSolver(verbose=True)
		set_bbox_min = np.array([-1, -1, -1], dtype=float)
		set_bbox_max = np.array([1, 1, 1], dtype=float)
		solve_options = {
			'bbox_min': set_bbox_min,
			'bbox_max': set_bbox_max,
		}
		phi = solver.compute_distance_to_mesh(V, F, options=solve_options)
		bbox_min, bbox_max = solver.get_bbox()
		assert np.linalg.norm(bbox_min - set_bbox_min) < 1e-8, 'bbox_min does not match!'
		assert np.linalg.norm(bbox_max - set_bbox_max) < 1e-8, 'bbox_max does not match!'
		assert (len(bbox_min) == 3) and (len(bbox_max) == 3), 'Grid corners should be 3D positions.'
		assert not np.any((bbox_max - bbox_min) < 0), (
			'Grid corners should be returned in order of minimal -> maximal node corner.'
		)

	def test_to_grid_array(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatGridSolver(verbose=True)
		solve_options = {}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		phi_grid = solver.to_grid_array(phi)
		assert len(phi_grid.shape) == 3, 'SDF on grid should be a 3D NumPy array.'

	def test_compute_distance_to_mesh(self) -> None:
		V, F = pp3d.read_polygon_mesh(os.path.join(asset_path, 'bunny_small.obj'))
		solver = shm.SignedHeatGridSolver(verbose=True)
		solve_options = {'resolution': np.array([32, 32, 32])}
		phi = solver.compute_distance_to_mesh(V=V, F=F, options=solve_options)
		assert len(phi.shape) == 1, 'SDF should be a 1D NumPy array.'
		# Make sure distance is close to naive distance.
		nx, ny, nz = solver.get_grid_resolution()
		bbox_min, bbox_max = solver.get_bbox()
		Q = grid_node_positions(nx, ny, nz, bbox_min, bbox_max)
		N = area_weighted_vertex_normals(V, F)
		signed_distances = approximate_signed_distance(Q, V, N)
		span = np.amax(signed_distances) - np.amin(signed_distances)
		residual = (phi - signed_distances) / span
		residual[np.isnan(residual)] = 0.0
		assert np.mean(residual) < 2e-2, 'SDF not close to approximate ground-truth.'

	def test_compute_distance_to_point_cloud(self) -> None:
		P, N = read_point_cloud(os.path.join(asset_path, 'bunny.pc'))
		solver = shm.SignedHeatGridSolver(verbose=True)
		solve_options = {'resolution': np.array([32, 32, 32])}
		phi = solver.compute_distance_to_point_cloud(P=P, N=N, options=solve_options)
		assert len(phi.shape) == 1, 'SDF should be a 1D NumPy array.'
		# Make sure distance is close to naive distance.
		nx, ny, nz = solver.get_grid_resolution()
		bbox_min, bbox_max = solver.get_bbox()
		Q = grid_node_positions(nx, ny, nz, bbox_min, bbox_max)
		signed_distances = approximate_signed_distance(Q, P, N)
		span = np.amax(signed_distances) - np.amin(signed_distances)
		residual = (phi - signed_distances) / span
		residual[np.isnan(residual)] = 0.0
		assert np.mean(residual) < 2e-2, 'SDF not close to approximate ground-truth.'
