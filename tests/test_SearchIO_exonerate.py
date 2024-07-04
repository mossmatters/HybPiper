"""
Tests for SearchIO ExonerateIO text parser.
"""

import os
import unittest
from Bio.SearchIO import read

TEST_DIR = "Exonerate"


class ExonerateTextCases(unittest.TestCase):
    fmt = "exonerate-text"

    def test_exonerate_results_frameshifts(self):
        """
        Test parsing exonerate output (exonerate_results_frameshifts.fasta).
        """

        exn_file = os.path.join(TEST_DIR, 'exonerate_results_frameshifts.fasta')
        qresult = read(exn_file, self.fmt)

        # Check common attributes:
        for hit in qresult:
            self.assertEqual(qresult.id, hit.query_id)
            for hsp in hit:
                self.assertEqual(hit.id, hsp.hit_id)
                self.assertEqual(qresult.id, hsp.query_id)

        self.assertEqual('Arecales_Arecaceae_Chamaedorea_seifrizii_Timilsena_et_al-4527', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual('exonerate', qresult.program)
        self.assertEqual('protein2genome:local', qresult.model)
        self.assertEqual(5, len(qresult))

        # Third hit
        hit = qresult[2]
        self.assertEqual('NODE_3_length_925_cov_11.824561', hit.id)
        self.assertEqual(1, len(hit.hsps))

        # First hit, only hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(4, len(hsp))
        self.assertEqual([], hsp.hit_split_codons)
        self.assertEqual(4, len(hsp.query_all))
        self.assertEqual(4, len(hsp.hit_all))
        self.assertEqual([1, 1, 1, 1], hsp.hit_strand_all)
        self.assertEqual([3, 1, 2, 3], hsp.hit_frame_all)  # i.e. has frameshifts
        self.assertEqual(123, hsp.query_start)
        self.assertEqual(302, hsp.query_end)
        self.assertEqual(203, hsp.hit_start)
        self.assertEqual(710, hsp.hit_end)
        self.assertEqual([(123, 158), (158, 209), (209, 273), (273, 302)], hsp.query_range_all)
        self.assertEqual([(203, 308), (309, 453), (454, 622), (623, 710)], hsp.hit_range_all)
        self.assertEqual([(158, 158), (209, 209), (273, 273)], hsp.query_inter_ranges)
        self.assertEqual([(308, 309), (453, 454), (622, 623)], hsp.hit_inter_ranges)
        self.assertEqual(
            'KYQTFTNPSDAKQYVKRQGAPIVVKADGLAAGKGV', hsp.query_all[0].seq[:40]
        )
        self.assertEqual(
            'QYQTFTNPSDAKKYIKEQGAPTVVKADGLAAGKGV', hsp.hit_all[0].seq[:40]
        )
        self.assertEqual(
            [':!!', '|||', '|||', '|||', '|||', '|||', '|||', '|||', '|||', '|||', '|||', '|||', ':!!', '|||',
             ':!!', '|||', '..!', '|||', '|||', '|||', '|||', '! !', '|||', '|||', '|||', '|||', '|||', '|||', '|||',
             '|||', '|||', '|||', '|||', '|||', '|||'],
            hsp[0].aln_annotation['similarity'][:40]
        )
        self.assertEqual(
            ['CAG', 'TAT', 'CAA', 'ACT', 'TTT', 'ACA', 'AAC', 'CCA', 'TCC', 'GAT', 'GCA', 'AAG', 'AAA', 'TAT',
             'ATT', 'AAG', 'GAG', 'CAA', 'GGA', 'GCG', 'CCC', 'ACT', 'GTT', 'GTC', 'AAG', 'GCT', 'GAT', 'GGG', 'TTG',
             'GCT', 'GCT', 'GGG', 'AAA', 'GGA', 'GTT'],
            hsp[0].aln_annotation['hit_annotation'][:40]
        )
