"""Port of tests/test_superposition.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka

ONE_MISMATCHING_EXPECTED = [
    (11.4030913638093, 11.4030913638093, 11.1658735912377),
    (12.1343637879364, 12.1343637879364, 12.5552911970792),
    (12.8656362120636, 12.8656362120636, 13.9447088029208),
    (13.5969086361907, 13.5969086361907, 15.3341264087623),
]

PLAIN_BACKBONE_EXPECTED = [
    (0.950860168366356, -4.10490242639581, 1.88352814333396),
    (0.50054990344796, -3.69993907500445, 0.514678133530478),
    (-0.347244679657512, -4.71592059085125, -0.00580232173625572),
    (-0.318861255763144, -2.43716579580966, 0.473326017357792),
    (0.506320258152502, -1.32662194302102, 0.113179074950623),
    (-0.0811260469759572, 0.154118620506023, 0.0681432878857751),
    (-0.647318283320744, 0.253840078930681, -1.40665838767353),
    (0.237737393186242, 0.145726060576478, -2.52768907793944),
    (-0.603024580420967, 0.282861734520597, -3.76182161161158),
    (-1.56591988716846, -0.783314033669751, -3.80479671565305),
    (-1.41697603269085, 1.56517005180424, -3.86385794226367),
    (-1.71799695715546, 1.88014731841389, -5.20722860018113),
]

EXTENDED_BACKBONE_EXPECTED = [
    (0.987110766185577, -4.09071771489897, 1.9236822102414),
    (0.498993067687067, -3.70426391350659, 0.562462326533447),
    (-0.388755973663699, -4.71109589650857, 0.0933514680094707),
    (-0.292926693772989, -2.42405654051827, 0.52155788233221),
    (0.54395097846502, -1.34065484814116, 0.109472290902603),
    (-1.17780200834469, -4.12249011848516, -0.933417051780487),
    (-2.43228362842333, -4.89306776930891, -0.85803171543564),
    (-3.01959843323915, -5.3413321995708, -1.99405078244382),
    (-0.011455206873104, 0.151981739114691, 0.052253520613716),
    (-0.625836124161192, 0.233795179516065, -1.40424936281698),
    (0.217366956113348, 0.0817923356419141, -2.55212814599476),
    (-0.662209557696309, 0.212387424607011, -3.75962973153933),
    (-1.64960487803313, -0.831943050574988, -3.7466166618369),
    (-1.45033183726805, 1.51072162058409, -3.86177956266596),
    (-1.79041850416729, 1.80432120889078, -5.2006631799918),
    (-2.93672182938542, -0.21171102176145, -3.72940666545558),
    (-3.83674604788897, -1.07392594379514, -2.92836744398792),
    (-4.64473104553373, -2.04674049128451, -3.42943939468367),
]


class TestSuperposition(LLKATestCase):
    def assertPoints(self, actual, expected, count: int | None = None) -> None:
        count = len(expected) if count is None else count
        for idx in range(count):
            self.assertPoint(actual[idx], expected[idx], f'point {idx}')

    def test_shifted(self) -> None:
        a = llka.Points([(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)])
        b = llka.Points([(11, 11, 11), (12, 12, 12), (13, 13, 13), (14, 14, 14)])

        rmsd = llka.superpose(a, b)

        self.assertPoints(a, list(b))
        self.assertFlt(rmsd, 0.0, 'RMSD is wrong')

    def test_shifted_2(self) -> None:
        a = llka.Points([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
        b = llka.Points([(10, 10, 10), (10, 11, 10), (10, 11, 11), (10, 10, 11)])

        rmsd = llka.superpose(a, b)

        self.assertPoints(a, list(b))
        self.assertFlt(rmsd, 0.0, 'RMSD is wrong')

    def test_one_mismatching(self) -> None:
        a = llka.Points([(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)])
        b = [(11, 11, 11), (12, 12, 12), (13, 13, 13), (14, 14, 17)]

        rmsd = llka.superpose(a, b)

        self.assertPoints(a, ONE_MISMATCHING_EXPECTED)
        self.assertFlt(rmsd, 1.0869242161333, 'RMSD is wrong')

    def test_mismatching_sizes(self) -> None:
        with self.assertRaises(llka.MismatchingSizesError):
            llka.superpose(llka.Points([(1, 1, 1), (2, 2, 2)]), [(1, 1, 1)])
        with self.assertRaises(llka.MismatchingSizesError):
            llka.rmsd([(1, 1, 1), (2, 2, 2)], [(1, 1, 1)])

    def test_real_plain_backbone(self) -> None:
        ref_bkbn = llka.extract_backbone(make_structure('REF_AB01_ATOMS'))
        real_bkbn = llka.extract_backbone(make_structure('REAL_1BNA_A_3_4_ATOMS'))

        rmsd = llka.superpose(real_bkbn, ref_bkbn)

        self.assertPoints([a.coords for a in real_bkbn], PLAIN_BACKBONE_EXPECTED)
        self.assertFlt(rmsd, 0.147657613575976, 'RMSD is wrong')

    def test_real_extended_backbone(self) -> None:
        ref_bkbn = llka.extract_extended_backbone(make_structure('REF_AB01_ATOMS'))
        real_bkbn = llka.extract_extended_backbone(make_structure('REAL_1BNA_A_3_4_ATOMS'))

        rmsd = llka.superpose(real_bkbn, ref_bkbn)

        self.assertPoints([a.coords for a in real_bkbn], EXTENDED_BACKBONE_EXPECTED, 14)
        self.assertFlt(rmsd, 0.176588293714412, 'RMSD is wrong')

    def test_real_extended_backbone_view(self) -> None:
        ref = make_structure('REF_AB01_ATOMS')
        ref_bkbn = llka.extract_extended_backbone_view(ref)
        real_bkbn = llka.extract_extended_backbone(make_structure('REAL_1BNA_A_3_4_ATOMS'))

        rmsd = llka.superpose(real_bkbn, ref_bkbn)

        self.assertPoints([a.coords for a in real_bkbn], EXTENDED_BACKBONE_EXPECTED, 14)
        self.assertFlt(rmsd, 0.176588293714412, 'RMSD is wrong')

    def test_transformation_matrix_points(self) -> None:
        a = llka.Points([(1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)])
        b = llka.Points([(11, 11, 11), (12, 12, 12), (13, 13, 13), (14, 14, 17)])

        transformation = llka.superposition_matrix(a, b)
        llka.apply_transformation(a, transformation)

        self.assertPoints(a, ONE_MISMATCHING_EXPECTED)

    def test_transformation_matrix_structure(self) -> None:
        ref_bkbn = llka.extract_extended_backbone(make_structure('REF_AB01_ATOMS'))
        real_bkbn = llka.extract_extended_backbone(make_structure('REAL_1BNA_A_3_4_ATOMS'))

        transformation = llka.superposition_matrix(real_bkbn, ref_bkbn)
        llka.apply_transformation(real_bkbn, transformation)

        self.assertPoints([a.coords for a in real_bkbn], EXTENDED_BACKBONE_EXPECTED, 14)

    def test_transformation_matrix_structure_view(self) -> None:
        ref = make_structure('REF_AB01_ATOMS')
        real = make_structure('REAL_1BNA_A_3_4_ATOMS')
        ref_bkbn = llka.extract_extended_backbone_view(ref)
        real_bkbn = llka.extract_extended_backbone_view(real)

        transformation = llka.superposition_matrix(real_bkbn, ref_bkbn)

        expected = [
            [0.6781837391735, -0.278187893249037, 0.680204610371256, -23.0397139352514],
            [0.362607977150515, -0.678396158822793, -0.638978956305519, 15.7850194094942],
            [0.639204404574515, 0.679992755655633, -0.359204094392586, -22.3197349277867],
            [0, 0, 0, 1],
        ]

        self.assertEqual(transformation.shape, (4, 4), 'Transformation matrix has invalid shape')
        for r in range(4):
            for c in range(4):
                self.assertFlt(transformation[r, c], expected[r][c], f'Wrong value in transformation matrix at [{r}; {c}]')

    def test_invalid_transformation_matrix(self) -> None:
        with self.assertRaises(llka.InvalidArgumentError):
            llka.apply_transformation(llka.Points([(1, 1, 1)]), llka.Matrix(3, 3))

    def test_centroid(self) -> None:
        pts = [(0, 0, 0), (2, 4, 6)]
        self.assertEqual(llka.centroid(pts), llka.Point(1, 2, 3))

        stru = make_structure('REF_AB01_ATOMS')
        expected = [sum(a.coords[i] for a in stru) / len(stru) for i in range(3)]
        self.assertPoint(llka.centroid(stru), expected)


if __name__ == '__main__':
    unittest.main()
