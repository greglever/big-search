#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the search module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import pytest
import logging

LOG = logging.getLogger(__name__)

try:
    GENOME_PATH = os.path.join(os.path.dirname(__file__), 'data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
    assert os.path.exists(GENOME_PATH)
except AssertionError:
    LOG.error('Human reference genome file not found. Download the human reference genome')
    raise


# @pytest.mark.parametrize(("pattern", 'mismatches', 'exphits'), [
#     (b'TGGATGTGAAATGAGTCAAG', 3, 'data/TGGATGTGAAATGAGTCAAG-results.sam'),
#     (b'GGGTGGGGGGAGTTTGCTCC', 3, 'data/vegfa-site1-results.sam'),
# ])
# def test_search(pattern, mismatches, exphits_path):
#     # TODO
#     result = set()
#     expected_hits = set()
#     with open(exphits_path, 'rb') as exphits:
#         for hit in exphits.readlines():
#             # TODO use pysam to parse the expected result records if needed.
#             expected_hits.add(hit)
#     # TODO implement a more details comparison function if needed
#     assert expected_hits.difference(result) is None


