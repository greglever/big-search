#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the search module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals


def alphabet_index(c, alphabet=None):
    """
    Returns the index of the given character in the English alphabet, counting from 0.
    """
    if alphabet is None:
        alphabet = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    # TODO: What if we encounter a part of the sequence that isn't A,C,G or T ?
    return alphabet[c]


def match_length(S, idx1, idx2):
    """
    Returns the length of the match of the substrings of S beginning at idx1 and idx2.
    """
    if idx1 == idx2:
        return len(S) - idx1
    match_count = 0
    while idx1 < len(S) and idx2 < len(S) and S[idx1] == S[idx2]:
        match_count += 1
        idx1 += 1
        idx2 += 1
    return match_count


def fundamental_preprocess(S):
    """
    Returns Z, the Fundamental Preprocessing of S. Z[i] is the length of the substring
    beginning at i which is also a prefix of S. This pre-processing is done in O(n) time,
    where n is the length of S.
    """
    if len(S) == 0: # Handles case of empty string
        return []
    if len(S) == 1: # Handles case of single-character string
        return [1]
    z = [0 for x in S]
    z[0] = len(S)
    z[1] = match_length(S, 0, 1)
    for i in range(2, 1+z[1]): # Optimization from exercise 1-5
        z[i] = z[1]-i+1
    # Defines lower and upper limits of z-box
    l = 0
    r = 0
    for i in range(2+z[1], len(S)):
        if i <= r: # i falls within existing z-box
            k = i-l
            b = z[k]
            a = r-i+1
            if b < a: # b ends within existing z-box
                z[i] = b
            else: # b ends at or after the end of the z-box, we need to do an explicit match to the right of the z-box
                z[i] = a+match_length(S, a, r+1)
                l = i
                r = i+z[i]-1
        else: # i does not reside within existing z-box
            z[i] = match_length(S, 0, i)
            if z[i] > 0:
                l = i
                r = i+z[i]-1
    return z


def get_bad_character_table(S):
    """
    Generates R for S, which is an array indexed by the position of some character c in the
    English alphabet. At that index in R is an array of length |S|+1, specifying for each
    index i in S (plus the index after S) the next location of character c encountered when
    traversing S from right to left starting at i. This is used for a constant-time lookup
    for the bad character rule in the Boyer-Moore string search algorithm, although it has
    a much larger size than non-constant-time solutions.
    """
    if len(S) == 0:
        return [[] for a in range(4)]
    R = [[-1] for a in range(4)]
    alpha = [-1 for a in range(4)]
    for i, c in enumerate(S):
        alpha[alphabet_index(c=c)] = i
        for j, a in enumerate(alpha):
            R[j].append(a)
    return R


def get_good_suffix_table(S):
    """
    Generates L for S, an array used in the implementation of the strong good suffix rule.
    L[i] = k, the largest position in S such that S[i:] (the suffix of S starting at i) matches
    a suffix of S[:k] (a substring in S ending at k). Used in Boyer-Moore, L gives an amount to
    shift P relative to T such that no instances of P in T are skipped and a suffix of P[:L[i]]
    matches the substring of T matched by a suffix of P in the previous match attempt.
    Specifically, if the mismatch took place at position i-1 in P, the shift magnitude is given
    by the equation len(P) - L[i]. In the case that L[i] = -1, the full shift table is used.
    Since only proper suffixes matter, L[0] = -1.
    """
    L = [-1 for c in S]
    N = fundamental_preprocess(S[::-1]) # S[::-1] reverses S
    N.reverse()
    for j in range(0, len(S)-1):
        i = len(S) - N[j]
        if i != len(S):
            L[i] = j
    return L


def get_full_shift_table(S):
    """
    Generates F for S, an array used in a special case of the good suffix rule in the Boyer-Moore
    string search algorithm. F[i] is the length of the longest suffix of S[i:] that is also a
    prefix of S. In the cases it is used, the shift magnitude of the pattern P relative to the
    text T is len(P) - F[i] for a mismatch occurring at i-1.
    """
    F = [0 for c in S]
    Z = fundamental_preprocess(S)
    longest = 0
    for i, zv in enumerate(reversed(Z)):
        longest = max(zv, longest) if zv == i+1 else longest
        F[-i-1] = longest
    return F


def string_search(pattern, text):
    """
    Implementation of the Boyer-Moore string search algorithm. This finds all occurrences of P
    in T, and incorporates numerous ways of pre-processing the pattern to determine the optimal
    amount to shift the string and skip comparisons. In practice it runs in O(m) (and even
    sublinear) time, where m is the length of T. This implementation performs a case-insensitive
    search on ASCII alphabetic characters, spaces not included.
    """
    if len(pattern) == 0 or len(text) == 0 or len(text) < len(pattern):
        return []

    matches = []

    # Preprocessing
    bad_character_table = get_bad_character_table(pattern)
    good_suffix_table = get_good_suffix_table(pattern)
    full_shift_table = get_full_shift_table(pattern)

    k = len(pattern) - 1      # Represents alignment of end of P relative to T
    previous_k = -1     # Represents alignment in previous phase (Galil's rule)
    while k < len(text):
        i = len(pattern) - 1  # Character to compare in P
        h = k           # Character to compare in T
        while i >= 0 and h > previous_k and pattern[i] == text[h]:   # Matches starting from end of P
            i -= 1
            h -= 1
        if i == -1 or h == previous_k:  # Match has been found (Galil's rule)
            matches.append(k - len(pattern) + 1)
            k += len(pattern) - full_shift_table[1] if len(pattern) > 1 else 1
        else:   # No match, shift by max of bad character and good suffix rules
            char_shift = i - bad_character_table[alphabet_index(text[h])][i]
            if i+1 == len(pattern):   # Mismatch happened on first attempt
                suffix_shift = 1
            elif good_suffix_table[i+1] == -1:   # Matched suffix does not appear anywhere in P
                suffix_shift = len(pattern) - full_shift_table[i + 1]
            else:               # Matched suffix appears in P
                suffix_shift = len(pattern) - good_suffix_table[i + 1]
            shift = max(char_shift, suffix_shift)
            previous_k = k if shift >= i+1 else previous_k  # Galil's rule
            k += shift
    return matches
