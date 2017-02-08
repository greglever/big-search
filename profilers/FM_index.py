
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions for creating FM index for text T
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

from itertools import combinations


def suffixArray(s):
    ''' Given T return suffix array SA(T).  Uses "sorted"
        function for simplicity, which is probably very slow. '''
    suffix_array = [(s[i:], i) for i in xrange(0, len(s))]
    """
    For an example Text ACGTTGCA the suffix array would be the following:
    [('ACGTTGCA', 0),
    ('CGTTGCA', 1),
    ('GTTGCA', 2),
    ('TTGCA', 3),
    ('TGCA', 4),
    ('GCA', 5),
    ('CA', 6),
    ('A', 7)]

    which for the human reference genome will give 4.5*10^18 rows...
    4,500,000,000,000,000,000.00

    so 1,500,000,000 times as big as the human genome

    np.save(file="./small_test_array", arr=np.ones([1, 3000000000]))
    creates a file that is 22GB...
    so 66,000,000,000 GB for all of them in total...
    """
    sorted_suffix_array = sorted(suffix_array)
    return map(lambda x: x[1], sorted_suffix_array)  # extract, return just offsets


def bwtFromSa(t, sa=None):
    ''' Given T, returns BWT(T) by way of the suffix array. '''
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si-1])
    return (''.join(bw), dollarRow) # return string-ized version of list bw


class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''

    def __init__(self, bw, cpIval=4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}  # checkpoints
        self.cpIval = cpIval  # spacing between checkpoints
        tally = {}  # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i in xrange(0, len(bw)):
            tally[bw[i]] += 1  # up to *and including*
            if (i % cpIval) == 0:
                for c in tally.iterkeys():
                    self.cps[c].append(tally[c])

    def rank(self, bw, c, row):
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc


class FmIndex:
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.

        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''

    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ''' Take only the suffix-array entries for every nth suffix.  Keep
            suffixes at offsets 0, n, 2n, etc with respect to the text.
            Return map from the rows to their suffix-array values. '''
        ssa = {}
        for i in xrange(0, len(sa)):
            # We could use i % n instead of sa[i] % n, but we lose the
            # constant-time guarantee for resolutions
            if sa[i] % n == 0:
                ssa[i] = sa[i]
        return ssa

    def __init__(self, t, cpIval=4, ssaIval=4):
        if t[-1] != '$':
            t += '$'  # add dollar if not there already
        #TODO(Greg): Make sure there is a $ at the end of the reference genome

        # Get BWT string and offset of $ within it
        sa = suffixArray(t)
        self.bwt, self.dollarRow = bwtFromSa(t, sa)


        # Get downsampled suffix array, taking every 1 out of 'ssaIval'
        # elements w/r/t T
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.slen = len(self.bwt)
        # Make rank checkpoints
        self.cps = FmCheckpoints(self.bwt, cpIval)
        # Calculate # occurrences of each character
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.iteritems()):
            self.first[c] = totc
            totc += count

    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.first:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.first.iterkeys()):
                if c < cc: return self.first[cc]
            return self.first[cc]
        else:
            return self.first[c]

    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1  # closed (inclusive) interval
        for i in xrange(len(p) - 1, -1, -1):  # from right to left
            l = self.cps.rank(self.bwt, p[i], l - 1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r) + self.count(p[i]) - 1
            if r < l:
                break
        return l, r + 1

    def resolve(self, row):
        ''' Given BWM row, return its offset w/r/t T '''

        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            return self.cps.rank(self.bwt, c, row - 1) + self.count(c)

        nsteps = 0
        while row not in self.ssa:
            row = stepLeft(row)
            nsteps += 1
        return self.ssa[row] + nsteps

    def hasSubstring(self, p):
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l

    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen - 1

    def occurrences(self, p):
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        return [self.resolve(x) for x in xrange(l, r)]

    def approximate_match(self, pattern, max_allowed_mismatches):
        for combination_number in xrange(len(pattern), len(pattern)-max_allowed_mismatches-1, -1):
            for pattern_combination in list(combinations(pattern, combination_number)):
                sub_pattern = ''.join([word for word in pattern_combination])
                print(sub_pattern, self.occurrences(p=sub_pattern))

    # TODO(Greg): Test this:
    #             In[3]: p, t = "CAT", "TTGTGTGCATGTTGTTTCATCATTTAGAGATACATTGCGCTGCATCATGGTCA"
    #
    #             In[4]: fm = FmIndex(t)
    #
    #             In[5]: fm.approximate_match(pattern=p, max_allowed_mismatches=2)
    #             (u'CAT', [42, 17, 45, 7, 32, 20])
    #             (u'CA', [51, 42, 17, 45, 7, 32, 20])
    #             (u'CT', [39])
    #             (u'AT', [29, 43, 18, 46, 8, 33, 21])
