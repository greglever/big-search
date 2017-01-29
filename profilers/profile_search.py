#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Profile the Matcher Module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

from .boyer_moore import string_search


def profile_genome_search():
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
    genome_file = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
    pattern = b'TGGATGTGAAATGAGTCAAG'
    testText = b'CGCTAAAAGCTAGAGCTACGCGACGATCAGCACTACGTGGATGTGAAATGAGTCAAGCGCGCTAGACGACTACGACTAGCAGCATCGATCGATCGATCG'
    maxmismatch = 3
    string_search(pattern=pattern, text=testText)
    # TODO: (Greg) Implementation


def main():
    profile_genome_search()

if __name__ == '__main__':
    main()
