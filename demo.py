import os, sys
import argparse
import numpy as np
import time

import polyscope as ps
import polyscope.imgui as psim

# Path to where the bindings live
sys.path.append(os.path.join(os.path.dirname(__file__), "build/"))
sys.path.append(os.path.join(os.path.dirname(__file__), "src/"))

import shm3d

def read_polygon_mesh(filepath):
	'''
	Read a surface mesh, assumed to be OBJ format.

	Args:
		filepath: string

	Returns:
		vertices: |V| x 3 NumPy array
		faces: list of lists; each sublist represents a face of arbitrary degree (with 0-indexed vertices)
	'''
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

def write_surface_mesh(vertices, faces, filepath):
	'''
	Write a surface mesh as an OBJ file.

	Args:
		vertices: an |V| x 3 NumPy array
		faces: a list of lists; each sublist represents a face of arbitrary degree (assumed to be 0-indexed)
		filepath: output filepath

	Returns:
		Nothing. Just writes to `filepath`.
	'''
	with open(filepath, 'w') as file:
		n_vertices = vertices.shape[0]
		n_faces = len(faces)
		for i in range(n_vertices):
			file.write('v %f %f %f\n' %(vertices[i,0], vertices[i,1], vertices[i,2]))
		for i in range(n_faces):
			f_idxs = ' '.join([str(v) + ' ' for v in faces[i]])
			file.write('f ' + f_idxs + '\n')

	return np.array(vertices, dtype=np.float64), faces

def read_point_cloud(filepath):
	'''
	Read a point cloud and its normals, assumed to be a plaintext file of newline-separated point positions and normals.

	Args:

	Returns: 
		positions: |P| x 3 NumPy array.
		normals: |P| x 3 NumPy array.
	'''
	points = []
	normals = []
	with open(filepath, "r") as file:
		for line in file:
			parts = line.strip().split()
			if not parts:
				continue
			elif parts[0] == "v":
				pos = list(map(float, parts[1:4]))
				points.append(pos)
			elif parts[0] == "vn":
				vec = list(map(float, parts[1:4]))
				normals.append(vec)

	return np.array(points), np.array(normals)

class DemoSolver():
	'''
	Handles general solving (from mesh or point cloud), as well as visualization.
	'''
	def __init__(self, input_mode, use_grid=False, h_coef=0, verbose=True, headless=False):

		self.input_mode = input_mode
		self.grid = use_grid
		self.h_coef = h_coef
		self.t_coef = 1.
		self.verbose = verbose
		self.headless = headless

		self.points = []
		self.point_normals = []
		self.vertices = []
		self.faces = []
		self.iso_vertices = []
		self.iso_faces = []

		self.rebuild = True
		self.constraint_mode = "ZeroSet"
		self.isoval = 0.
		self.phi = []
		self.contoured = False
		self.last_solver_mode = "tet"
		self.output_dir = "output/"
		self.ps_plane = ""
		self.grid_scalar_q = []

		self.tet_solver = SignedHeatTetSolver(verbose)
		self.grid_solver = SignedHeatGridSolver(verbose)

	def contour(self):

		if self.last_solver_mode == "tet":
			self.iso_vertices, self.iso_faces = self.tet_solver.isosurface(self.phi, self.isoval)
			ps.register_surface_mesh("isosurface", self.iso_vertices, self.iso_faces)
		else:
			self.grid_scalar_q.set_isosurface_level(self.isoval)
			self.grid_scalar_q.set_isosurface_viz_enabled(True)
			self.grid_scalar_q.set_slice_planes_affect_isosurface(False)
			self.grid_scalar_q.register_isosurface_as_mesh("isosurface")

		self.contoured = True;
		ps.get_surface_mesh("isosurface").set_ignore_slice_plane(self.ps_plane, True);

	def callback(self):

		if psim.Button("Solve"):
			solve()
		if psim.RadioButton("on tet mesh", self.mesh_mode == "tet"): 
			self.mesh_mode == "tet"
		if psim.RadioButton("on grid", self.mesh_mode == "grid"): 
			self.mesh_mode == "grid"

		_, self.t_coef = psim.InputFloat("tCoef (diffusion time)", self.t_coef)
		self.rebuild, self.h_coef = psim.InputFloat("hCoef (mesh spacing)", self.h_coef)

		if self.mesh_mode != "grid":
			if psim.RadioButton("Constrain zero set", self.constraint_mode == "ZeroSet"): 
				self.constraint_mode = "ZeroSet"
			if psim.RadioButton("Constrain multiple levelsets", self.constraint_mode == "Multiple"): 
				self.constraint_mode = "Multiple"
			if psim.RadioButton("No level set constraints", self.constraint_mode == "None"): 
				self.constraint_mode = "None"

		if len(self.phi) > 0:
			psim.Separator()
			psim.Text("Contour options")
			psim.Separator()
			if psim.SliderFloat("Contour (drag slider)", self.isoval, np.min(self.phi), np.max(self.phi)): 
				contour()
			if psim.InputFloat("Contour (enter value)", self.isoval): 
				contour()

			if self.contoured:
				if psim.Button("Export isosurface"):
					if self.last_solver_mode == "grid":
						psIsoMesh = ps.get_surface_mesh("isosurface")
						isoFilename = self.output_dir + "/isosurface_" + str(self.isoval) + ".obj"
						write_surface_mesh(positions, psIsoMesh.vertices, psIsoMesh.faces)

	def solve(self):

		cmap = "viridis"
		if self.mesh_mode == "tet":
			if self.verbose: print("\nSolving on tet mesh...")
			t1 = time.time()
			self.phi = self.tet_solver.compute_distance_to_mesh(self.vertices, self.faces, self.constraint_mode, self.t_coef, self.h_coef, self.rebuild, self.scale) if (self.input_mode == "mesh") else self.tet_solver.compute_distance_to_point_cloud(self.points, self.normals, self.constraint_mode, self.t_coef, self.h_coef, self.rebuild, self.scale)
			t2 = time.time()
			if self.verbose: print("Solve time (s): %f" %(t2-t1))
			if not self.headless:
				ps.get_volume_mesh("domain").add_vertex_scalar_quantity("GSD", self.phi).set_color_map(cmap).set_isolines_enabled(True).set_enabled(True)
		else:
			if self.verbose: print("\nSolving on grid...")
			t1 = time.time()
			self.phi = self.grid_solver.compute_distance_to_mesh(self.vertices, self.faces, self.t_coef, self.h_coef, self.rebuild, self.scale) if (self.input_mode == "mesh") else self.grid_solver.compute_distance_to_point_cloud(self.points, self.normals, self.t_coef, self.h_coef, self.rebuild, self.scale)
			t2 = time.time()
			if self.verbose: print("Solve time (s): %f" %(t2-t1))
			if not self.headless:
				self.grid_scalar_q = ps.get_volume_mesh("domain").add_node.scalar_quantity("GSD", self.phi).set_color_map(cmap).set_isolines_enabled(True).setEnabled(True)

		if self.verbose:
			print("min: %f \tmax: %f" %(np.min(self.phi), np.max(self.phi)))

		if not headless:
			ps.remove_last_scene_slice_plane()
			self.ps_plane = ps.add_scene_slice_plane()
			self.ps_plane.set_draw_plane(False)
			self.ps_plane.set_draw_widget(True)
			if self.mesh_mode == "tet":
				self.ps_plane.set_volume_mesh_to_inspect("domain")
			if self.input_mode == "mesh":
				ps.get_surface_mesh("mesh").set_ignore_slice_plane(self.ps_plane, True)
			else:
				ps.get_point_cloud("point cloud").set_ignore_slice_plane(self.ps_plane, True)

		self.last_solver_mode = self.mesh_mode
		self.rebuild = False


def main():

	parser = argparse.ArgumentParser("signed-heat")
	parser.add_argument("input", help="A mesh or point cloud file.", type=str)
	parser.add_argument('-h') # tet/grid spacing
	parser.add_argument('-g', '--grid', action='store_true')
	parser.add_argument('-v', '--verbose', action='store_true')
	parser.add_argument('-l', '--headless', action='store_true')
	args = parser.parse_args()

	if not args.input:
		raise RuntimeError("Please specify an input curve or mesh.")

	filepath = os.path.basename(args.input)
	meshname, ext = os.path.splitext(filepath)
	input_mode = "cloud" if (ext == ".pc") else "mesh" 
	demo_solver = DemoSolver(input_mode, args.grid, args.h, args.verbose, args.headless)
	if ext != ".pc":
		demo_solver.vertices, demo_solver.faces = read_polygon_mesh(filepath)
	else:
		demo_solver.points, demo_solver.point_normals = read_point_cloud(filepath)


	if not args.headless:
		ps.init()
		ps.set_user_callback(demo_solver.callback)
		if ext != ".pc":
			ps.register_surface_mesh("mesh", demo_solver.vertices, demo_solver.faces)
		else:
			ps.register_point_cloud("point cloud", demo_solver.points)

		ps.show()
	else:
		demo_solver.solve()

if __name__ == '__main__':
	main()