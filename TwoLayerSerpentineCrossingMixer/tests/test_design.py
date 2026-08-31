from __future__ import annotations

import importlib.util
import copy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SPEC = importlib.util.spec_from_file_location("multifidelity_design", ROOT / "multifidelity_design.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

CAD_SPEC = importlib.util.spec_from_file_location(
    "m10_cad", ROOT / "FlowCase" / "two_layer_serpentine_crossing_cad.py"
)
CAD_MODULE = importlib.util.module_from_spec(CAD_SPEC)
assert CAD_SPEC.loader is not None
CAD_SPEC.loader.exec_module(CAD_MODULE)

PILOT_SPEC = importlib.util.spec_from_file_location(
    "multifidelity_pilot", ROOT / "multifidelity_pilot.py"
)
PILOT_MODULE = importlib.util.module_from_spec(PILOT_SPEC)
assert PILOT_SPEC.loader is not None
PILOT_SPEC.loader.exec_module(PILOT_MODULE)


class MultifidelityDesignTests(unittest.TestCase):
    def test_interlayer_aperture_manifest(self) -> None:
        config = copy.deepcopy(MODULE.load_yaml(MODULE.REFERENCE_PATH))
        config["geometry"]["number_of_units"] = 2
        _, manifest = CAD_MODULE.build_geometry(config)
        connections = manifest["interlayer_connections"]
        expected_vertical_area = 3 * 0.15e-3 * 1.07e-3
        self.assertAlmostEqual(
            connections["vertical_segment_open_area_m2"], expected_vertical_area
        )
        self.assertGreater(
            connections["crossing_open_area_outside_vertical_segments_m2"], 0.0
        )
        self.assertAlmostEqual(
            connections["total_open_area_m2"],
            connections["vertical_segment_open_area_m2"]
            + connections["crossing_open_area_outside_vertical_segments_m2"],
        )

    def test_campaign_contract(self) -> None:
        campaign = MODULE.load_yaml(MODULE.CAMPAIGN_PATH)
        MODULE.validate_campaign(campaign)
        self.assertEqual(
            campaign["scientific_framing"]["classification"],
            "m10_inspired_reconstruction",
        )
        self.assertFalse(campaign["scientific_framing"]["exact_reproduction_claim"])
        self.assertEqual(float(campaign["design"]["operating_reynolds_number"]), 10.0)
        self.assertEqual(float(campaign["fidelities"]["coarse"]["nominal_cell_size_m"]), 2.4e-5)
        self.assertEqual(float(campaign["fidelities"]["fine"]["nominal_cell_size_m"]), 1.3e-5)

    def test_pilot_is_deterministic_subset_of_initialization(self) -> None:
        campaign = PILOT_MODULE.campaign_config()
        first_pool = PILOT_MODULE.sobol_pool(campaign)
        second_pool = PILOT_MODULE.sobol_pool(campaign)
        self.assertEqual(first_pool, second_pool)
        anchors = PILOT_MODULE.paired_anchor_indices(campaign, first_pool)
        self.assertEqual(len(first_pool), 24)
        self.assertEqual(len(anchors), 6)
        self.assertEqual(len(set(anchors)), 6)
        self.assertTrue(all(0 <= index < len(first_pool) for index in anchors))

    def test_requested_scalar_schemes_are_configured(self) -> None:
        research = MODULE.load_yaml(ROOT / "research_config.yaml")
        numerics = research["scalar_transport_numerics"]
        self.assertEqual(numerics["transported_field"], "T")
        self.assertEqual(numerics["convection_scheme"], "Gauss linearUpwind gradT")
        self.assertEqual(
            numerics["gradient_scheme"], "cellLimited pointCellsLeastSquares 1"
        )
        schemes = (ROOT / "ScalarTransportCase" / "system" / "fvSchemes").read_text()
        self.assertIn("div(phi,T)      Gauss linearUpwind gradT;", schemes)
        self.assertIn("gradT           cellLimited pointCellsLeastSquares 1;", schemes)

    def test_reference_geometry_ratios(self) -> None:
        params = MODULE.reference_parameters()
        self.assertAlmostEqual(params["H_over_P"], 1.07 / 0.64)
        self.assertAlmostEqual(params["w_over_P"], 0.30 / 0.64)
        self.assertAlmostEqual(params["D_over_P"], 0.30 / 0.64)
        self.assertAlmostEqual(params["b_over_P"], 0.15 / 0.64)

    def test_same_design_has_distinct_fidelity_inputs(self) -> None:
        params = MODULE.reference_parameters()
        coarse = MODULE.materialize_design(params, "coarse")
        fine = MODULE.materialize_design(params, "fine")
        self.assertEqual(coarse["geometry"], fine["geometry"])
        self.assertEqual(coarse["fidelity"]["coordinate"], 0.0)
        self.assertEqual(fine["fidelity"]["coordinate"], 1.0)
        self.assertGreater(fine["fidelity"]["nominal_relative_cost"], coarse["fidelity"]["nominal_relative_cost"])
        self.assertLess(fine["fidelity"]["cell_size_m"], coarse["fidelity"]["cell_size_m"])

    def test_lower_bound_design_is_mesh_resolved(self) -> None:
        campaign = MODULE.load_yaml(MODULE.CAMPAIGN_PATH)
        lower = {
            name: float(bounds["lower"])
            for name, bounds in campaign["design"]["parameters"].items()
        }
        MODULE.materialize_design(lower, "coarse")
        MODULE.materialize_design(lower, "fine")


if __name__ == "__main__":
    unittest.main()
