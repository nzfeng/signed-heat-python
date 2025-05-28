#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "signed_heat_grid_solver.h"
#include "signed_heat_tet_solver.h"

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace geometrycentral::pointcloud;

SignedHeat3DOptions toSignedHeatOptions(const std::string& levelSetConstraint, double tCoef, double hCoef, bool rebuild,
                                        double scale) {

  auto toLower = [&](const std::string& s) -> std::string {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return std::tolower(c); });
    return t;
  };

  SignedHeat3DOptions options;
  if (toLower(levelSetConstraint) == "none") {
    options.levelSetConstraint = LevelSetConstraint::None;
  }
  if (toLower(levelSetConstraint) == "zeroset") {
    options.levelSetConstraint = LevelSetConstraint::ZeroSet;
  }
  if (toLower(levelSetConstraint) == "multiple") {
    options.levelSetConstraint = LevelSetConstraint::Multiple;
  }
  options.tCoef = tCoef;
  options.hCoef = hCoef;
  options.rebuild = rebuild;
  options.scale = scale;
  return options;
}

std::unique_ptr<VertexPositionGeometry> makeSurfaceGeometry(const DenseMatrix<double>& vertices,
                                                            const std::vector<std::vector<size_t>>& faces) {
  std::vector<Vector3> vertexPositions;
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = makeSurfaceMeshAndGeometry(faces, vertexPositions);
  return std::move(geometry);
}

std::unique_ptr<PointPositionNormalGeometry> makePointCloudGeometry(const DenseMatrix<double>& positions,
                                                                    const DenseMatrix<double>& normals) {

  size_t nPts = positions.rows();
  PointCloud cloud = PointCloud(nPts);
  PointData<Vector3> pointPositions = PointData<Vector3>(cloud);
  PointData<Vector3> pointNormals = PointData<Vector3>(cloud);
  for (size_t i = 0; i < nPts; i++) {
    for (int j = 0; j < 3; j++) {
      pointPositions[i][j] = positions(i, j);
      pointNormals[i][j] = normals(i, j);
    }
  }
  std::unique_ptr<PointPositionNormalGeometry> pointGeom = std::unique_ptr<PointPositionNormalGeometry>(
      new PointPositionNormalGeometry(cloud, pointPositions, pointNormals));
  return std::move(pointGeom);
}

// A wrapper class for SignedHeatTetSolver, which exposes the parameters of `options`, and passes mesh data as Eigen
// arrays.
class SignedHeatTetSolverWrapper {

public:
  SignedHeatTetSolverWrapper(bool verbose) {
    solver.reset(new SignedHeatTetSolver());
    solver->VERBOSE = verbose;
  }

  Vector<double> compute_distance_to_mesh(const DenseMatrix<double>& vertices,
                                          const std::vector<std::vector<size_t>>& faces, std::string levelSetConstraint,
                                          double tCoef, double hCoef, bool rebuild, double scale) {

    std::unique_ptr<VertexPositionGeometry> geometry = makeSurfaceGeometry(vertices, faces);
    SignedHeat3DOptions options = toSignedHeatOptions(levelSetConstraint, tCoef, hCoef, rebuild, scale);
    return solver->computeDistance(*geometry, options);
  }

  Vector<double> compute_distance_to_point_cloud(const DenseMatrix<double>& points, const DenseMatrix<double>& normals,
                                                 std::string levelSetConstraint, double tCoef, double hCoef,
                                                 bool rebuild, double scale) {

    std::unique_ptr<PointPositionNormalGeometry> pointGeom = makePointCloudGeometry(points, normals);
    SignedHeat3DOptions options = toSignedHeatOptions(levelSetConstraint, tCoef, hCoef, rebuild, scale);
    return solver->computeDistance(*pointGeom, options);
  }

  std::tuple<DenseMatrix<double>, std::vector<std::vector<size_t>>> isosurface(const Vector<double>& phi,
                                                                               const double& isoval) {
    std::unique_ptr<SurfaceMesh> isoMesh;
    std::unique_ptr<VertexPositionGeometry> isoGeom;
    solver->isosurface(isoMesh, isoGeom, phi, isoval);
    DenseMatrix<double> vertices(isoMesh->nVertices(), 3);
    for (size_t i = 0; i < isoMesh->nVertices(); i++) {
      for (int j = 0; j < 3; j++) {
        vertices(i, j) = isoGeom->vertexPositions[i][j];
      }
    }
    return std::make_tuple(vertices, isoMesh->getFaceVertexList());
  }

private:
  std::unique_ptr<SignedHeatTetSolver> solver;
};

// A wrapper class for SignedHeatGridSolver, which exposes the parameters of `options`, and passes mesh data as Eigen
// arrays.
class SignedHeatGridSolverWrapper {

public:
  SignedHeatGridSolverWrapper(bool verbose) {
    solver.reset(new SignedHeatGridSolver());
    solver->VERBOSE = verbose;
  }

  Vector<double> compute_distance_to_mesh(const DenseMatrix<double>& vertices,
                                          const std::vector<std::vector<size_t>>& faces, double tCoef, double hCoef,
                                          bool rebuild, double scale) {

    std::unique_ptr<VertexPositionGeometry> geometry = makeSurfaceGeometry(vertices, faces);
    SignedHeat3DOptions options = toSignedHeatOptions("None", tCoef, hCoef, rebuild, scale);
    return solver->computeDistance(*geometry, options);
  }

  Vector<double> compute_distance_to_point_cloud(const DenseMatrix<double>& points, const DenseMatrix<double>& normals,
                                                 double tCoef, double hCoef, bool rebuild, double scale) {

    std::unique_ptr<PointPositionNormalGeometry> pointGeom = makePointCloudGeometry(points, normals);
    SignedHeat3DOptions options = toSignedHeatOptions("None", tCoef, hCoef, rebuild, scale);
    return solver->computeDistance(*pointGeom, options);
  }

private:
  std::unique_ptr<SignedHeatGridSolver> solver;
};


// binding code
// clang-format off

NB_MODULE(shm3d_bindings, m) {

	nb::class_<SignedHeatTetSolverWrapper>(m, "SignedHeatTetSolver")
      .def(nb::init<bool>())
      .def("compute_distance_to_mesh", &SignedHeatTetSolverWrapper::compute_distance_to_mesh, 
      	nb::arg("vertices"), 
      	nb::arg("faces"),
      	nb::arg("level_set_constraint"),
      	nb::arg("t_coef"),
      	nb::arg("h_coef"),
      	nb::arg("rebuild"),
      	nb::arg("scale"))
      .def("compute_distance_to_point_cloud", &SignedHeatTetSolverWrapper::compute_distance_to_point_cloud, 
      	nb::arg("points"),
      	nb::arg("normals"),
        nb::arg("level_set_constraint"),
      	nb::arg("t_coef"),
      	nb::arg("h_coef"),
      	nb::arg("rebuild"),
      	nb::arg("scale"))
      .def("isosurface", &SignedHeatTetSolverWrapper::isosurface,
        nb::arg("phi"),
        nb::arg("isoval"));

  nb::class_<SignedHeatGridSolverWrapper>(m, "SignedHeatGridSolver")
      .def(nb::init<bool>())
      .def("compute_distance_to_mesh", &SignedHeatGridSolverWrapper::compute_distance_to_mesh, 
      	nb::arg("vertices"), 
      	nb::arg("faces"),
      	nb::arg("t_coef"),
      	nb::arg("h_coef"),
      	nb::arg("rebuild"),
      	nb::arg("scale"))
      .def("compute_distance_to_point_cloud", &SignedHeatGridSolverWrapper::compute_distance_to_point_cloud, 
      	nb::arg("points"),
      	nb::arg("normals"),
      	nb::arg("t_coef"),
      	nb::arg("h_coef"),
      	nb::arg("rebuild"),
      	nb::arg("scale"));
}