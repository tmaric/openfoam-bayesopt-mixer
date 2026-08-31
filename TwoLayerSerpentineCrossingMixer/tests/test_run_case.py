from __future__ import annotations

import importlib.util
from pathlib import Path
import csv
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("run_case", ROOT / "run_case.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class RunCaseTests(unittest.TestCase):
    def test_square_duct_pressure_reference(self) -> None:
        mu = 8.8e-4
        length = 3.55e-3
        side = 3.0e-4
        velocity = 0.044
        pressure = MODULE.fully_developed_rectangular_duct_pressure_drop_pa(
            mu, length, side, side, velocity
        )
        # Darcy Poiseuille number f*Re for a square duct is about 56.91.
        poiseuille = 2.0 * pressure * side**2 / (mu * velocity * length)
        self.assertAlmostEqual(poiseuille, 56.91, places=2)

    def test_review_lead_split_default_and_override(self) -> None:
        source = ROOT / "FlowCase" / "two_layer_serpentine_crossing_cad.yaml"
        with tempfile.TemporaryDirectory() as directory:
            default_path = Path(directory) / "default.yaml"
            default = MODULE.prepare_geometry_config("review", source, default_path)
            self.assertEqual(default["geometry"]["number_of_units"], 6)
            self.assertAlmostEqual(default["geometry"]["inlet_lead_length"], 0.08)
            self.assertAlmostEqual(default["geometry"]["outlet_lead_length"], 0.08)

            audit_path = Path(directory) / "audit.yaml"
            audit = MODULE.prepare_geometry_config(
                "review", source, audit_path, inlet_lead_mm=0.02, outlet_lead_mm=0.14
            )
            self.assertAlmostEqual(audit["geometry"]["inlet_lead_length"], 0.02)
            self.assertAlmostEqual(audit["geometry"]["outlet_lead_length"], 0.14)
            self.assertEqual(yaml.safe_load(audit_path.read_text()), audit)

    def test_snappy_background_is_materialized_without_absolute_paths(self) -> None:
        manifest = {
            "patches": {
                "inlet1": {"bounds_m": {"xmin": 0.0, "xmax": 0.0, "ymin": 0.0, "ymax": 1.0e-4, "zmin": 0.0, "zmax": 1.0e-4}},
                "inlet2": {"bounds_m": {"xmin": 0.0, "xmax": 0.0, "ymin": -1.0e-4, "ymax": 0.0, "zmin": -1.0e-4, "zmax": 0.0}},
                "outlet": {"bounds_m": {"xmin": 1.0e-3, "xmax": 1.0e-3, "ymin": -1.0e-4, "ymax": 1.0e-4, "zmin": -1.0e-4, "zmax": 1.0e-4}},
                "walls": {"bounds_m": {"xmin": 0.0, "xmax": 1.0e-3, "ymin": -2.0e-4, "ymax": 2.0e-4, "zmin": -1.0e-4, "zmax": 1.0e-4}},
            },
            "derived": {"snappy_location_in_mesh_m": [5.0e-5, 5.0e-5, 5.0e-5]},
        }
        with tempfile.TemporaryDirectory() as directory:
            flow = Path(directory)
            system = flow / "system"
            system.mkdir()
            for name in ("blockMeshDict", "snappyHexMeshDict"):
                source = ROOT / "FlowCase" / "system" / name
                (system / name).write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
            report = MODULE.materialize_snappy_background(
                flow, manifest, 2.5e-5, meshing_ranks=4
            )
            materialized = "\n".join(
                (system / name).read_text(encoding="utf-8")
                for name in ("blockMeshDict", "snappyHexMeshDict")
            )
        self.assertNotIn("__", materialized)
        self.assertNotIn("/home/", materialized)
        self.assertEqual(report["generator"], "snappyHexMesh")
        self.assertEqual(report["meshing_ranks"], 4)

    def test_boundary_gate_rejects_nonempty_background_patch(self) -> None:
        boundary = """
3
(
inlet
{ type patch; nFaces 10; startFace 100; }
outlet
{ type patch; nFaces 10; startFace 110; }
walls
{ type wall; nFaces 40; startFace 120; }
background
{ type patch; nFaces 1; startFace 160; }
)
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "boundary"
            path.write_text(boundary, encoding="utf-8")
            with self.assertRaises(RuntimeError):
                MODULE.validate_boundary_topology(path)

    def test_mesh_gate_accepts_hex_dominant_fine_cut_cells_without_tets(self) -> None:
        log = """
    faces:            1480000
    cells:            486100
    hexahedra:        468000
    prisms:           17600
    tet wedges:       1
    tetrahedra:       0
    polyhedra:        499
   *There are 8 faces with concave angles between consecutive edges. Max concave angle = 21.608069 degrees.
 ***Concave cells (using face planes) found, number of cells: 8
Checking faces in error :
    non-orthogonality > 55 degrees : 0
    pyramid volume < 1e-30 : 0
    tet quality < 1e-15 : 0
    concavity > 30 degrees : 0
    skewness > 2 : 0
    interpolation weights < 0.05 : 0
    volume ratio < 0.02 : 0
    face twist < 0.02 : 0
    determinant < 0.01 : 0

Failed 1 mesh checks.
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log.checkMesh"
            path.write_text(log, encoding="utf-8")
            report = MODULE.validate_mesh_quality(path, "fine")
        self.assertGreater(report["hexahedral_cell_fraction"], 0.95)
        self.assertEqual(report["tetrahedra"], 0)
        self.assertLess(report["concave_cell_fraction"], 2.0e-5)

    def test_mesh_gate_rejects_polyhedral_failure(self) -> None:
        log = """
    faces:            1480000
    cells:            486100
    hexahedra:        480000
    prisms:           5800
    tet wedges:       0
    tetrahedra:       0
    polyhedra:        300
 ***Concave cells (using face planes) found, number of cells: 3203
Checking faces in error :
    non-orthogonality > 55 degrees : 2
    pyramid volume < 1e-30 : 0
    tet quality < 1e-15 : 0
    concavity > 30 degrees : 18
    skewness > 2 : 4
    interpolation weights < 0.05 : 0
    volume ratio < 0.02 : 0
    face twist < 0.02 : 0
    determinant < 0.01 : 0

Failed 2 mesh checks.
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log.checkMesh"
            path.write_text(log, encoding="utf-8")
            with self.assertRaises(RuntimeError):
                MODULE.validate_mesh_quality(path, "fine")

    def test_scalar_history_gate(self) -> None:
        fieldnames = [
            "flux_weighted_intensity_of_segregation",
            "flux_weighted_mean_concentration",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mixing.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fieldnames)
                writer.writeheader()
                for index in range(60):
                    writer.writerow(
                        {
                            fieldnames[0]: 0.1 + 1.0e-7 * index,
                            fieldnames[1]: 0.5 + 1.0e-8 * index,
                        }
                    )
            report = MODULE.validate_scalar_history(path, window=50)
        self.assertLess(report["intensity_span"], 1.0e-4)

    def test_scalar_history_gate_rejects_relative_drift(self) -> None:
        fieldnames = [
            "flux_weighted_intensity_of_segregation",
            "flux_weighted_mean_concentration",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mixing.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fieldnames)
                writer.writeheader()
                for index in range(60):
                    writer.writerow(
                        {
                            fieldnames[0]: 1.0e-4 + 1.0e-6 * index,
                            fieldnames[1]: 0.5,
                        }
                    )
            with self.assertRaises(RuntimeError):
                MODULE.validate_scalar_history(path, window=50)

    def test_scalar_bounds_gate(self) -> None:
        text = "min(T) = -0.0002 in cell 1 max(T) = 1.0003 in cell 2\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log.scalarBounds"
            path.write_text(text, encoding="utf-8")
            report = MODULE.validate_scalar_bounds(path)
        self.assertAlmostEqual(report["undershoot"], 2.0e-4)
        self.assertAlmostEqual(report["overshoot"], 3.0e-4)

    def test_scalar_bounds_gate_rejects_excursion_above_declared_cap(self) -> None:
        text = "min(T) = -0.0021 in cell 1 max(T) = 1.0 in cell 2\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "log.scalarBounds"
            path.write_text(text, encoding="utf-8")
            with self.assertRaises(RuntimeError):
                MODULE.validate_scalar_bounds(path)

    def test_protocols_keep_distinct_diffusivities(self) -> None:
        config = MODULE.load_yaml(MODULE.RESEARCH_CONFIG)
        self.assertEqual(config["software"]["openfoam_version"], "v2606")
        self.assertEqual(config["software"]["mesh_generator"], "snappyHexMesh")
        original = config["literature_reproduction"]["fluid"]["scalar_diffusivity_m2_s"]
        review = config["review_matched_benchmark"]["scalar_diffusivity_m2_s"]
        self.assertEqual(float(original), 1.0e-11)
        self.assertEqual(float(review), 1.0e-10)


if __name__ == "__main__":
    unittest.main()
