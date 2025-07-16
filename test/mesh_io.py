import numpy as np


def read_polygon_mesh(filepath: str) -> tuple[np.ndarray, list[list[int]]]:
	"""
	Read a surface mesh, assumed to be OBJ format.

	Args:
		filepath: string

	Returns:
		vertices: |V| x 3 NumPy array
		faces: list of lists; each sublist represents a face of arbitrary degree (with 0-indexed vertices)
	"""
	vertices = []
	faces = []
	with open(filepath, 'r') as file:
		for line in file:
			parts = line.strip().split()
			if not parts:
				continue
			elif parts[0] == 'v':
				vertex = list(map(float, parts[1:4]))
				vertices.append(vertex)
			elif parts[0] == 'f':
				face = [int(p.split('/')[0]) - 1 for p in parts[1:]]
				faces.append(face)

	return np.array(vertices, dtype=np.float64), faces


def write_surface_mesh(vertices: np.ndarray, faces: list[list[int]], filepath: str) -> None:
	"""
	Write a surface mesh as an OBJ file.

	Args:
		vertices: an |V| x 3 NumPy array
		faces: a list of lists; each sublist represents a face of arbitrary degree (assumed to be 0-indexed)
		filepath: output filepath

	Returns:
		Nothing. Just writes to `filepath`.
	"""
	with open(filepath, 'w') as file:
		n_vertices = vertices.shape[0]
		n_faces = len(faces)
		for i in range(n_vertices):
			file.write('v %f %f %f\n' % (vertices[i, 0], vertices[i, 1], vertices[i, 2]))
		for i in range(n_faces):
			f_idxs = ' '.join([str(v + 1) + ' ' for v in faces[i]])  # OBJs are 1-indexed
			file.write('f ' + f_idxs + '\n')

	return np.array(vertices, dtype=np.float64), faces


def read_point_cloud(filepath: str) -> tuple[np.ndarray, np.ndarray]:
	"""
	Read a point cloud and its normals, assumed to be a plaintext file of newline-separated point positions and normals.

	Args:
		filepath: string

	Returns:
		positions: |P| x 3 NumPy array.
		normals: |P| x 3 NumPy array.
	"""
	points = []
	normals = []
	with open(filepath, 'r') as file:
		for line in file:
			parts = line.strip().split()
			if not parts:
				continue
			elif parts[0] == 'v':
				pos = list(map(float, parts[1:4]))
				points.append(pos)
			elif parts[0] == 'vn':
				vec = list(map(float, parts[1:4]))
				normals.append(vec)

	return np.array(points), np.array(normals)
