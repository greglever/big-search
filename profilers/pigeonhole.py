#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

from .boyer_moore import string_search


def approximate_match(pattern, text, max_allowed_mismatches):
    pattern_length = len(pattern)
    text_length = len(text)
    segment_length = int(round(pattern_length / (max_allowed_mismatches + 1)))
    all_matches = set()
    for segment_index in range(max_allowed_mismatches + 1):
        start = segment_index + segment_length
        finish = min((segment_index + 1)*segment_length, (pattern_length))
        match_positions = string_search(pattern=pattern[start:finish], text=text)
        for match_position in match_positions:
            if match_position < start or match_position - start + pattern_length > text_length:
                continue  # fall out of the loop if this evaluates as true

            mismatch_count = 0
            for index_value in range(start):
                if pattern[index_value] != text[match_position - start + index_value]:
                    mismatch_count += 1  # we've encountered a mismatch
                    if mismatch_count > max_allowed_mismatches:
                        break  # This partition of the pattern won't yield appromximate matches within allowance
            for suffix_index in range(finish, pattern_length):
                if pattern[suffix_index] != text[match_position - start + suffix_index]:
                    mismatch_count += 1
                    if mismatch_count > max_allowed_mismatches:
                        break

            if mismatch_count <= max_allowed_mismatches:
                all_matches.add(match_position - start)
    return list(all_matches)
