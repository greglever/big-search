#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Profile the Matcher Module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import sys

from profilers.pigeonhole import approximate_match


def profile_genome_search(pattern, max_allowed_mismatches):
    """
    Profile the sort bucket step
    
    The idea is to use the "pidgeonhole principle" which will allow us
    to  use the Boyer-Moore algorithm to find approximate (as opposed to exact) matches by 
    dividing the pattern P up into partitions, looking for each of 
    those partitions using an exact matching algorithm (Boyer-Moore). And then, 
    wherever we find an occurrence of one of those partitions, 
    we can do a verification step where we look in the neighborhood to 
    see if we see a full occurrence of P with up to the maximum number 
    of differences allowed. 
    """
    # genome_file = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
    testText = 'CGCTAAAAGCTAGAGCTACGCGACGATCAGCACTACGTGGATGTGAAATGAGTCAAGCGCGCTAGACGACTACGACTAGCAGCATCGATTGGATGGGAAAGGAGTCAAGCGATCGATCG'
    print(approximate_match(pattern=pattern, text=testText, max_allowed_mismatches=max_allowed_mismatches))
    # TODO: (Greg) Implementation


if __name__ == '__main__':
    # pattern = 'TGGATGTGAAATGAGTCAAG'
    if len(sys.argv) < 3:
        print("Usage: profile_search.py [pattern] [maximum allowed number of mismatches]")
        print("You must supply pattern and maximium allowed mismatch arguments to profile_search.py ")
        exit()
    pattern = sys.argv[1]
    max_allowed_mismatches = int(sys.argv[2])
    profile_genome_search(pattern=pattern, max_allowed_mismatches=max_allowed_mismatches)
