"""Port of tests/test_classification.cpp."""

from __future__ import annotations

import math
import unittest

from _support import ASSETS_DIR, LLKATestCase, make_structure

import llka
from llka import ClassificationViolation as V


def make_context() -> llka.ClassificationContext:
    def load(name: str, kind: llka.ResourceType) -> llka.ResourceList:
        return llka.load_resource_file(ASSETS_DIR / name, kind)

    limits = llka.ClassificationLimits(
        minimum_nearest_neighbors=7,
        number_of_used_nearest_neighbors=11,
        average_neighbors_torsion_cutoff=math.radians(28.0),
        nearest_neighbor_torsions_cutoff=math.radians(28.0),
        total_distance_cutoff=math.radians(60.0),
        pseudorotation_cutoff=math.radians(72.0),
        minimum_cluster_votes=0.001111,
    )

    return llka.ClassificationContext(
        clusters=load('clusters.csv', llka.ResourceType.CLUSTERS),
        golden_steps=load('golden_steps.csv', llka.ResourceType.GOLDEN_STEPS),
        confals=load('confals.csv', llka.ResourceType.CONFALS),
        cluster_nu_angles=load('nu_angles.csv', llka.ResourceType.AVERAGE_NU_ANGLES),
        confal_percentiles=load('confal_percentiles.csv', llka.ResourceType.CONFAL_PERCENTILES),
        limits=limits,
        max_close_enough_rmsd=0.5,
    )


class TestClassification(LLKATestCase):
    ctx: llka.ClassificationContext

    @classmethod
    def setUpClass(cls) -> None:
        cls.ctx = make_context()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.ctx.close()

    def test_sugar_pucker_naming(self) -> None:
        self.assertEqual(llka.name_to_sugar_pucker('C1end'), llka.SugarPucker.C1_ENDO, 'wrong sugar pucker')
        self.assertEqual(llka.name_to_sugar_pucker("C3'exo"), llka.SugarPucker.C3_EXO, 'wrong sugar pucker')
        self.assertEqual(llka.name_to_sugar_pucker("C4' exo"), llka.SugarPucker.C4_EXO, 'wrong sugar pucker')
        self.assertEqual(llka.name_to_sugar_pucker("O4' endo"), llka.SugarPucker.O4_ENDO, 'wrong sugar pucker')
        self.assertEqual(llka.name_to_sugar_pucker("O1'exo"), llka.SugarPucker.O4_EXO, 'wrong sugar pucker')
        with self.assertRaises(ValueError):
            llka.name_to_sugar_pucker('C3 exo')

        self.assertEqual(llka.sugar_pucker_to_name(llka.SugarPucker.C1_ENDO, llka.SugarPuckerNameBrevity.TERSE), 'C1endo')
        self.assertEqual(llka.sugar_pucker_to_name(llka.SugarPucker.C1_EXO, llka.SugarPuckerNameBrevity.FANCY), "C1' exo")
        self.assertEqual(llka.sugar_pucker_to_name(llka.SugarPucker.INVALID, llka.SugarPuckerNameBrevity.VERY_TERSE), '')

    def test_classify_simple(self) -> None:
        step = self.ctx.classify_step(make_structure('REAL_1BNA_A_1_2_ATOMS'))

        self.assertEqual(step.violations, V.OK, 'detected classification violations where there were not supposed to be any')
        self.assertTrue(step.is_assigned)
        self.assertEqual(step.assigned_ntc, llka.NtC.BB04, 'wrong NtC class')
        self.assertEqual(step.closest_golden_step, '1mjo_F_DG_10_DA_11', 'wrong golden step name')
        self.assertFlt(step.rmsd_to_closest_ntc, 0.366977877615573, 'RMSD to closest NtC exceeds tolerance')

        self.assertFlt(step.confal_score.CC, 80.5998065953948, 'wrong CC confal score')
        self.assertFlt(step.confal_score.NN, 99.9952502902733, 'wrong NN confal score')
        self.assertFlt(step.confal_score.total, 32.9484362098594, 'wrong total confal score')

    def test_classify_thorough(self) -> None:
        step = self.ctx.classify_step(make_structure('REAL_3VOK_U_1_2_ATOMS'))

        self.assertEqual(step.violations, V.OK, 'detected classification violations where there were not supposed to be any')
        self.assertEqual(step.assigned_ntc, llka.NtC.BBS1, 'wrong NtC class')
        self.assertEqual(step.assigned_cana, llka.CANA.SYN, 'wrong CANA class')
        self.assertEqual(step.assigned_ntc, step.closest_ntc, 'assigned and closest NtCs do not match')
        self.assertEqual(step.assigned_cana, step.closest_cana, 'assigned and closest CANAs do not match')
        self.assertEqual(step.closest_golden_step, '1lv5_F_DA_1_DC_2', 'wrong golden step')

        self.assertFlt(step.confal_score.CC, 35.6080290838014, 'wrong CC confal score')
        self.assertFlt(step.confal_score.NN, 43.1082983683617, 'wrong NN confal score')
        self.assertFlt(step.confal_score.total, 72.3729761097013, 'wrong total confal score')

        self.assertEqual(step.sugar_pucker_1, llka.SugarPucker.C2_ENDO, 'wrong sugar pucker on first nucleotide')
        self.assertEqual(step.sugar_pucker_2, llka.SugarPucker.C1_EXO, 'wrong sugar pucker on second nucleotide')
        self.assertFlt(step.ribose_pseudorotation_1, 2.84035778887486, 'wrong ribose pseudorotation on first ribose')
        self.assertFlt(step.tau_2, 0.7674388643588, 'wrong tau_2')
        self.assertFlt(step.nu_angles_2.nu_2, -0.608524908531341, 'wrong nu2 angle on second nucleotide')

        avg = self.ctx.average_confal([step])
        self.assertFlt(avg.score, 72.372976109701284, 'wrong average confal score')
        self.assertFlt(avg.percentile, 94.599999999999994, 'wrong confal percentile')

    def test_classify_violations(self) -> None:
        step = self.ctx.classify_step(make_structure('REAL_1DK1_B_5_6_ATOMS'))

        self.assertEqual(
            step.violations,
            V.NOT_ENOUGH_NEAREST_NEIGHBORS
            | V.AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT
            | V.NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT
            | V.BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES,
            'wrong violations',
        )
        self.assertFalse(step.is_assigned)
        self.assertIsNone(step.assigned_ntc, 'no NtC should be assigned')
        self.assertEqual(step.closest_ntc, llka.NtC.IC03, 'wrong NtC class')

    def test_classify_not_classifiable(self) -> None:
        with self.assertRaises(llka.MultipleAltIdsError):
            self.ctx.classify_step(make_structure('REAL_1DK1_B_27_28_ATOMS'))

    def test_classify_multiple(self) -> None:
        results = self.ctx.classify_steps([
            make_structure('REAL_1BNA_A_1_2_ATOMS'),
            make_structure('REAL_1DK1_B_27_28_ATOMS'),
            make_structure('REAL_3VOK_U_1_2_ATOMS'),
        ])

        self.assertEqual(len(results), 3)
        self.assertTrue(results[0].ok)
        self.assertEqual(results[0].step.assigned_ntc, llka.NtC.BB04)
        self.assertFalse(results[1].ok)
        self.assertIsNone(results[1].step)
        self.assertIsInstance(results[1].error, llka.MultipleAltIdsError)
        self.assertEqual(results[2].step.assigned_ntc, llka.NtC.BBS1)

        avg = self.ctx.average_confal(results)
        expected = (results[0].step.confal_score.total + results[2].step.confal_score.total) / 2
        self.assertFlt(avg.score, expected, 'wrong average confal score')

        with self.assertRaises(llka.NothingToClassifyError):
            self.ctx.classify_steps([])

    def test_get_cluster(self) -> None:
        cluster = self.ctx.cluster_for_ntc(llka.NtC.AA00)
        self.assertEqual(cluster.ntc, llka.NtC.AA00, 'wrong cluster NtC')

        with self.assertRaises(llka.InvalidArgumentError):
            self.ctx.cluster_for_ntc(llka.NtC.INVALID)

    def test_get_confal(self) -> None:
        cluster = self.ctx.cluster_for_ntc(llka.NtC.BB00)
        confal = self.ctx.confal_for_ntc(llka.NtC.BB00)
        self.assertEqual(confal.cluster_number, cluster.number)

    def test_confal_percentile(self) -> None:
        self.assertFlt(self.ctx.confal_percentile(72.372976109701284), 94.599999999999994)
        with self.assertRaises(ValueError):
            self.ctx.confal_percentile(101.0)
        with self.assertRaises(ValueError):
            self.ctx.confal_percentile(math.nan)

    def test_violation_names(self) -> None:
        self.assertEqual(
            llka.classification_violation_to_name(V.CC_TOO_LOW),
            'LLKA_CLASSIFICATION_E_CC_TOO_LOW',
        )

    def test_context_from_python_objects(self) -> None:
        def load(name: str, kind: llka.ResourceType) -> list:
            return list(llka.load_resource_file(ASSETS_DIR / name, kind))

        limits = llka.ClassificationLimits(7, 11, math.radians(28.0), math.radians(28.0), math.radians(60.0),
                                           math.radians(72.0), 0.001111)
        with llka.ClassificationContext(
            clusters=load('clusters.csv', llka.ResourceType.CLUSTERS),
            golden_steps=load('golden_steps.csv', llka.ResourceType.GOLDEN_STEPS),
            confals=load('confals.csv', llka.ResourceType.CONFALS),
            cluster_nu_angles=load('nu_angles.csv', llka.ResourceType.AVERAGE_NU_ANGLES),
            confal_percentiles=load('confal_percentiles.csv', llka.ResourceType.CONFAL_PERCENTILES),
            limits=limits,
            max_close_enough_rmsd=0.5,
        ) as ctx:
            step = ctx.classify_step(make_structure('REAL_1BNA_A_1_2_ATOMS'))
            self.assertEqual(step.assigned_ntc, llka.NtC.BB04)

        with self.assertRaises(ValueError):
            ctx.classify_step(make_structure('REAL_1BNA_A_1_2_ATOMS'))

    def test_bad_limits(self) -> None:
        limits = llka.ClassificationLimits(7, 11, 0.0, 0.0, 0.0, 0.0, 0.0)
        with self.assertRaises(llka.BadClassificationLimitsError):
            llka.ClassificationContext(
                clusters=llka.load_resource_file(ASSETS_DIR / 'clusters.csv', llka.ResourceType.CLUSTERS),
                golden_steps=llka.load_resource_file(ASSETS_DIR / 'golden_steps.csv', llka.ResourceType.GOLDEN_STEPS),
                confals=llka.load_resource_file(ASSETS_DIR / 'confals.csv', llka.ResourceType.CONFALS),
                cluster_nu_angles=llka.load_resource_file(ASSETS_DIR / 'nu_angles.csv', llka.ResourceType.AVERAGE_NU_ANGLES),
                confal_percentiles=llka.load_resource_file(
                    ASSETS_DIR / 'confal_percentiles.csv', llka.ResourceType.CONFAL_PERCENTILES
                ),
                limits=limits,
                max_close_enough_rmsd=0.5,
            )

    def test_missing_resource_file(self) -> None:
        with self.assertRaises(llka.NoFileError):
            llka.load_resource_file(ASSETS_DIR / 'does_not_exist.csv', llka.ResourceType.CLUSTERS)
        with self.assertRaises(FileNotFoundError):
            llka.load_resource_file(ASSETS_DIR / 'does_not_exist.csv', llka.ResourceType.CLUSTERS)


if __name__ == '__main__':
    unittest.main()
