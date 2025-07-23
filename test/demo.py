import os
import argparse
import numpy as np
import time
import potpourri3d as pp3d

import polyscope as ps
import polyscope.imgui as psim

from mesh_io import *

import signed_heat_method as shm


class DemoSolver:
	"""
	Handles general solving (from mesh or point cloud), as well as visualization.
	"""

	def __init__(
		self,
		input_mode: str,
		use_grid: bool = False,
		bbox_min: np.ndarray = np.array([1.0, 1.0, 1.0]),
		bbox_max: np.ndarray = np.array([-1.0, -1.0, -1.0]),
		resolution: np.ndarray = np.array([32, 32, 32]),
		verbose: bool = True,
		headless: bool = False,
	) -> None:
		self.input_mode = input_mode
		self.mesh_mode = 'grid' if use_grid else 'tet'

		self.solve_options = {
			'level_set_constraint': 'ZeroSet',
			'rebuild': True,
			't_coef': 1.0,
			'bbox_min': bbox_min,
			'bbox_max': bbox_max,
			'resolution': resolution,
			'verbose': verbose,
		}
		self.headless = headless

		self.points = None
		self.point_normals = None
		self.vertices = None
		self.faces = None
		self.iso_vertices = None
		self.iso_faces = None

		self.isoval = 0.0
		self.phi = None
		self.contoured = False
		self.last_solver_mode = 'tet'
		self.output_dir = 'output/'
		self.ps_plane = None
		self.grid_scalar_q = None

		self.tet_solver = shm.SignedHeatTetSolver(verbose)
		self.grid_solver = shm.SignedHeatGridSolver(verbose)

	def contour(self) -> None:
		if self.last_solver_mode == 'tet':
			self.iso_vertices, self.iso_faces = self.tet_solver.isosurface(self.phi, self.isoval)
			ps.register_surface_mesh('isosurface', self.iso_vertices, self.iso_faces, enabled=True)
		else:
			self.grid_scalar_q.set_isosurface_level(self.isoval)
			self.grid_scalar_q.set_isosurface_viz_enabled(True)
			self.grid_scalar_q.set_slice_planes_affect_isosurface(False)
			self.grid_scalar_q.register_isosurface_as_mesh('isosurface')

		self.contoured = True
		ps.get_surface_mesh('isosurface').set_ignore_slice_plane(self.ps_plane, True)

	def solve(self) -> None:
		cmap = 'viridis'
		if self.mesh_mode != self.last_solver_mode:
			self.solve_options['rebuild'] = True
		if self.mesh_mode == 'tet':
			if self.solve_options['verbose']:
				print('\nSolving on tet mesh...')
			t1 = time.time()
			if self.input_mode == 'mesh':
				self.phi = self.tet_solver.compute_distance_to_mesh(
					V=self.vertices, F=self.faces, options=self.solve_options
				)
			else:
				self.phi = self.tet_solver.compute_distance_to_point_cloud(
					P=self.points, N=self.point_normals, options=self.solve_options
				)
			t2 = time.time()
			if self.solve_options['verbose']:
				print('Solve time (s): %f' % (t2 - t1))
			if not self.headless:
				if self.solve_options['rebuild']:
					ps.register_volume_mesh(
						'tet domain',
						self.tet_solver.get_vertices(),
						self.tet_solver.get_tets(),
					)
				# TODO: isolines not yet bound in Polyscope: https://github.com/nmwsharp/polyscope-py/issues/36
				ps.get_volume_mesh('tet domain').add_scalar_quantity('GSD', self.phi, cmap=cmap, enabled=True)
		else:
			if self.solve_options['verbose']:
				print('\nSolving on grid...')
			t1 = time.time()
			if self.input_mode == 'mesh':
				self.phi = self.grid_solver.compute_distance_to_mesh(
					V=self.vertices, F=self.faces, options=self.solve_options
				)
			else:
				self.phi = self.grid_solver.compute_distance_to_point_cloud(
					P=self.points, N=self.point_normals, options=self.solve_options
				)
			t2 = time.time()
			if self.solve_options['verbose']:
				print('Solve time (s): %f' % (t2 - t1))
			if not self.headless:
				if self.solve_options['rebuild']:
					grid_sizes = self.grid_solver.get_grid_resolution()
					bboxMin, bboxMax = self.grid_solver.get_bbox()
					ps.register_volume_grid('grid domain', grid_sizes, bboxMin, bboxMax)
				self.grid_scalar_q = ps.get_volume_grid('grid domain').add_scalar_quantity(
					'GSD',
					self.grid_solver.to_grid_array(self.phi),
					cmap=cmap,
					isolines_enabled=True,
					enabled=True,
				)

		if self.solve_options['verbose']:
			print('min: %f \tmax: %f' % (np.min(self.phi), np.max(self.phi)))

		if not self.headless:
			ps.remove_last_scene_slice_plane()
			self.ps_plane = ps.add_scene_slice_plane()
			self.ps_plane.set_draw_plane(False)
			self.ps_plane.set_draw_widget(True)
			if self.mesh_mode == 'tet':
				self.ps_plane.set_volume_mesh_to_inspect('tet domain')
			if self.input_mode == 'mesh':
				ps.get_surface_mesh('mesh').set_ignore_slice_plane(self.ps_plane, True)
			else:
				ps.get_point_cloud('point cloud').set_ignore_slice_plane(self.ps_plane, True)

		self.last_solver_mode = self.mesh_mode
		self.solve_options['rebuild'] = False

	def callback(self) -> None:
		if psim.Button('Solve'):
			self.solve()
		if psim.RadioButton('on tet mesh', self.mesh_mode == 'tet'):
			self.mesh_mode = 'tet'
		if psim.RadioButton('on grid', self.mesh_mode == 'grid'):
			self.mesh_mode = 'grid'

		_, self.solve_options['t_coef'] = psim.InputFloat('tCoef (diffusion time)', self.solve_options['t_coef'])
		changed_x, self.solve_options['resolution'][0] = psim.InputFloat(
			'Resolution (x-axis)', self.solve_options['resolution'][0]
		)
		if changed_x:
			self.solve_options['rebuild'] = True
		changed_y, self.solve_options['resolution'][1] = psim.InputFloat(
			'Resolution (y-axis)', self.solve_options['resolution'][1]
		)
		if changed_y:
			self.solve_options['rebuild'] = True
		changed_z, self.solve_options['resolution'][2] = psim.InputFloat(
			'Resolution (z-axis)', self.solve_options['resolution'][2]
		)
		if changed_z:
			self.solve_options['rebuild'] = True

		if self.mesh_mode != 'grid':
			if psim.RadioButton(
				'Constrain zero set',
				self.solve_options['level_set_constraint'] == 'ZeroSet',
			):
				self.solve_options['level_set_constraint'] = 'ZeroSet'
			if psim.RadioButton(
				'Constrain multiple levelsets',
				self.solve_options['level_set_constraint'] == 'Multiple',
			):
				self.solve_options['level_set_constraint'] = 'Multiple'
			if psim.RadioButton(
				'No level set constraints',
				self.solve_options['level_set_constraint'] == 'None',
			):
				self.solve_options['level_set_constraint'] = 'None'

		# if len(self.phi) > 0:
		# 	psim.Separator()
		# 	psim.Text("Contour options")
		# 	psim.Separator()
		# 	if psim.SliderFloat("Contour (drag slider)", self.isoval, np.min(self.phi), np.max(self.phi)):
		# 		self.contour()
		# 	if psim.InputFloat("Contour (enter value)", self.isoval):
		# 		self.contour()

		# 	if self.contoured:
		# 		if psim.Button("Export isosurface"):
		# 			if self.last_solver_mode == "grid":
		# 				psIsoMesh = ps.get_surface_mesh("isosurface")
		# 				isoFilename = self.output_dir + "/isosurface_" + str(self.isoval) + ".obj"
		# 				write_surface_mesh(psIsoMesh.vertices, psIsoMesh.faces, isoFilename)


def main() -> None:
	parser = argparse.ArgumentParser('signed-heat')
	parser.add_argument('input', help='A mesh or point cloud file.', type=str)
	parser.add_argument('-s', default=np.array([32, 32, 32], dtype=np.int64))  # tet/grid spacing
	parser.add_argument('-b', default=np.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0], dtype=np.float64))  # bbox parameters
	parser.add_argument('-g', '--grid', action='store_true', default=False)
	parser.add_argument('-v', '--verbose', action='store_true', default=True)
	parser.add_argument('-l', '--headless', action='store_true', default=False)
	args = parser.parse_args()

	if not args.input:
		raise RuntimeError('Please specify an input mesh or point cloud.')

	filename = os.path.basename(args.input)
	meshname, ext = os.path.splitext(filename)
	input_mode = 'cloud' if (ext == '.pc') else 'mesh'
	demo_solver = DemoSolver(
		input_mode=input_mode,
		use_grid=args.grid,
		bbox_min=args.b[:3],
		bbox_max=args.b[3:],
		resolution=args.s,
		verbose=args.verbose,
		headless=args.headless,
	)

	if ext != '.pc':
		demo_solver.vertices, demo_solver.faces = pp3d.read_polygon_mesh(args.input)
	else:
		demo_solver.points, demo_solver.point_normals = read_point_cloud(args.input)

	if not args.headless:
		ps.init()
		ps.set_user_callback(demo_solver.callback)
		if ext != '.pc':
			ps.register_surface_mesh('mesh', demo_solver.vertices, demo_solver.faces)
		else:
			ps.register_point_cloud('point cloud', demo_solver.points)

		ps.show()
	else:
		demo_solver.solve()


if __name__ == '__main__':
	main()
