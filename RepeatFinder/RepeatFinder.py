#!/usr/bin/pypy

######################################################################################
# RepeatFinder.py                                                                    #
#                                                                                    #
# SHELL : ./RepeatFinder.py <text file> <string type DNA/RNA> <minimum repeat size>  #
# PYTHON: from RepeatFinder import RepeatFinder                                      #
#                                                                                    #
# A Python module that outputs all repeats in a list of DNA/RNA sequences. When it's #
# imported, returns a dictionary of repeats (two types supported). Can also be used  #
# as a commandline PyPy application.                                                 #
#                                                                                    #
# Author: Ayaan Hossain (ain.hoss07@gmail.com)                                       #
######################################################################################

from string      import maketrans
from collections import defaultdict, OrderedDict, Counter
from itertools   import imap, islice, izip, count
import sys

class KMP(object):

    def __init__(self):
        self.text       = None
        self.pattern    = None
        self.start      = None
        self.transition = None

    def _set_transition_table(self):
        self.transition = [0] * len(self.pattern)
        for suffix_index in xrange(1, len(self.pattern)):
            prefix_index = self.transition[suffix_index - 1]
            while prefix_index > 0 and self.pattern[prefix_index] != self.pattern[suffix_index]:
                prefix_index = self.transition[prefix_index - 1]
            prefix_index += 1 if self.pattern[prefix_index] == self.pattern[suffix_index] else 0
            self.transition[suffix_index] = prefix_index

    def _get_exact_match(self):
        text_index, pattern_index = self.start, 0
        while text_index < len(self.text):
            if self.text[text_index] == self.pattern[pattern_index]:
                text_index, pattern_index = text_index + 1, pattern_index + 1
                if pattern_index == len(self.pattern):
                    yield text_index - pattern_index
                    pattern_index = self.transition[pattern_index - 1]
            else:
                if pattern_index == 0:
                    text_index += 1
                else:
                    pattern_index = self.transition[pattern_index - 1]

    def match_iter(self, text, pattern, start=0):
        self.text    = text
        self.pattern = pattern
        self.start   = start
        self._set_transition_table()
        return self._get_exact_match()

    def match(self, text, pattern, start=0):
        return list(self.match_iter(text, pattern, start))

class RepeatFinder(object):

    def __init__(self, str_type='DNA'):
        self.matcher     = KMP()
        self.compl_table = maketrans('ATGC', 'TACG') if str_type == 'DNA' else maketrans('AUGC', 'UACG')

    def _get_rev_comp(self, seq):
        return seq.translate(self.compl_table)[::-1]

    def _stream_seq_kmer_indices(self, seq, k):
        return (i for i in xrange(len(seq) - k + 1))

    def _is_left_maximal(self, seq_list, seq_k_id, k_start, k, j_start_list, rep_type):
        k_end = k_start + k
        if rep_type == 'direct' and k_start == 0:
            return True
        elif rep_type == 'invert' and k_end == len(seq_list[seq_k_id]):
            return True
        for seq_j_id, j_start_entry in j_start_list:
            if j_start_entry[0] == 0:
                return True
        for seq_j_id, j_start_entry in j_start_list:
            for j_start in j_start_entry:
                if rep_type == 'direct':
                    if seq_list[seq_k_id][k_start - 1] != seq_list[seq_j_id][j_start - 1]:
                        return True
                elif rep_type == 'invert':
                    if self._get_rev_comp(seq_list[seq_k_id][k_end]) != seq_list[seq_j_id][j_start - 1]:
                        return True
        return False

    def _is_right_maximal(self, seq_list, seq_k_id, k_start, k, j_start_list, rep_type):
        k_end = k_start + k
        if rep_type == 'direct' and k_end == len(seq_list[seq_k_id]):
            return True
        elif rep_type == 'invert' and k_start == 0:
            return True
        for seq_j_id, j_start_entry in j_start_list:
            if j_start_entry[-1] + k == len(seq_list[seq_j_id]):
                return True
        for seq_j_id, j_start_entry in j_start_list:
            for j_start in j_start_entry:
                if rep_type == 'direct':
                    if seq_list[seq_k_id][k_end] != seq_list[seq_j_id][j_start + k]:
                        return True
                elif rep_type == 'invert':
                    if self._get_rev_comp(seq_list[seq_k_id][k_start - 1]) != seq_list[seq_j_id][j_start + k]:
                        return True
        return False

    def _update_covered_kmers(self, covered_kmers, kmer):
        covered_kmers.add(kmer)

    def _get_maximal_repeat_indices(self, repeat_dict, covered_kmers, seq_list, seq_k_id, k_start, k, kmer, rep_type):
        if not kmer in covered_kmers:
            if kmer in repeat_dict:
                if rep_type in repeat_dict[kmer]:
                    return None
            j_start_list = []
            for seq_j_id in xrange(seq_k_id, len(seq_list)):
                if seq_k_id == seq_j_id:
                    j_start_entry = self.matcher.match(text=seq_list[seq_k_id], pattern=kmer, start=k_start+1)
                else:
                    j_start_entry = self.matcher.match(text=seq_list[seq_j_id], pattern=kmer, start=0)
                if j_start_entry:
                    j_start_list.append((seq_j_id, j_start_entry))
            if j_start_list:
                self._update_covered_kmers(covered_kmers, kmer)
                if self._is_left_maximal(seq_list, seq_k_id, k_start, k, j_start_list, rep_type):
                    if self._is_right_maximal(seq_list, seq_k_id, k_start, k, j_start_list, rep_type):
                        return j_start_list
        return None

    def _stream_all_maximal_repeats(self, repeat_dict, seq_k_id, seq_list, k_low, k_high):
        covered_kmers = True
        for k in xrange(k_low, k_high):
            covered_kmers = set()
            for k_start in self._stream_seq_kmer_indices(seq_list[seq_k_id], k):
                _kmer  = seq_list[seq_k_id][k_start:k_start + k]
                _rkmer = self._get_rev_comp(seq=_kmer)
                for kmer, rep_type in [(_kmer,'direct'), (_rkmer, 'invert')]:
                    kmer_repeat_indices = self._get_maximal_repeat_indices(repeat_dict, covered_kmers, seq_list, seq_k_id, k_start, k, kmer, rep_type)
                    if kmer_repeat_indices:
                        yield kmer, seq_k_id, [k_start], rep_type
                        for seq_j_id, j_start_entry in kmer_repeat_indices:
                            yield kmer, seq_j_id, j_start_entry, 'direct'
            if not covered_kmers:
                break

    def get_repeat_dict(self, seq_list, k_low, verbose=False, k_high=None):
        repeat_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
        pr_count = 0
        if verbose:
            print '{:>6}\t{:>15}\t{:>7}\t{:>7}\t  {:<9}'.format(' Length', 'Repeat', 'Type', 'Seq_ID', 'Locations')
        for seq_k_id in xrange(len(seq_list)):
            current_k_high = len(seq_list[seq_k_id]) - 1 if k_high is None else k_high
            prev_kmer = ''
            for kmer, seq_id, locs, rep_type in self._stream_all_maximal_repeats(repeat_dict, seq_k_id, seq_list, k_low, current_k_high):
                if verbose:
                    if kmer == prev_kmer:
                        print_kmer = len_kmer = print_rep = ''
                    else:
                        prev_kmer  = kmer
                        print_kmer = '...'.join([kmer[:5], kmer[-5:]]) if len(kmer) > 10 else kmer
                        len_kmer   = len(kmer)
                        print_rep  = rep_type
                        print ' --------------------------------------------------'
                        pr_count += 1
                    print ' {:>6}\t{:>15}\t{:>7}\t{:>7}\t  {:<9}'.format(len_kmer, print_kmer, print_rep, seq_id, ', '.join(imap(str, sorted(locs))))
                repeat_dict[kmer][rep_type][seq_id].update(locs)
        return repeat_dict

    def get_alt_repeat_dict(self, seq_list, k_low, verbose=False, k_high=None):
        repeat_dict = self.get_repeat_dict(seq_list, k_low, verbose, k_high)
        alt_dict    = defaultdict(lambda: defaultdict(set))
        for kmer in repeat_dict:
            for rep_type in repeat_dict[kmer]:
                for seq_id in repeat_dict[kmer][rep_type]:
                    alt_dict[kmer][seq_id].update(izip(repeat_dict[kmer][rep_type][seq_id], rep_type*len(repeat_dict[kmer][rep_type][seq_id])))
        repeat_dict.clear()
        return alt_dict

def main():
    if len(sys.argv) < 4:
        print ' Usage: ./RepeatFinder.py <text file> <string type DNA/RNA> <minimum repeat size>'
    else:
        with open(sys.argv[1]) as in_file:
            seq_list   = in_file.read().strip().split('\n')
        str_type       = sys.argv[2]
        k_low          = int(sys.argv[3])
        rep_finder_obj = RepeatFinder(str_type)
        rep_finder_obj.get_repeat_dict(seq_list, k_low, verbose=True)

if __name__ == '__main__':
    main()
