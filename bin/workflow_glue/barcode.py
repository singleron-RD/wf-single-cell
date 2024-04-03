#!/usr/bin/env python

import argparse
import itertools
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
import logging


def get_seq_str(seq, sub_pattern):
    """
    join seq slices.

    Args:
        seq: usually R1 read
        sub_pattern: [slice(0,8),slice(16,24)]

    Returns:
        joined intervals seq

    Raises:
        IndexError: if seq length is not enough.

    >>> sub_pattern_dict = [slice(0,2)]
    >>> seq = "A" * 2 + "T" * 2
    >>> get_seq_str(seq, sub_pattern_dict)
    'AA'
    """
    seq_len = len(seq)
    expect_len = sub_pattern[-1].stop
    if seq_len < expect_len:
        raise IndexError(f"read length({seq_len} bp) less than expected length({expect_len} bp) in read: {seq}")
    return "".join([seq[x] for x in sub_pattern])


def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch set.

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in itertools.combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in itertools.product(*seq_locs):
            seq_set.add(''.join(poss))
    return seq_set

def get_mismatch_dict(seq_list, n_mismatch=1):
    """
    Return:
    mismatch dict. Key: mismatch seq, value: seq in seq_list

    >>> seq_list = ["AACGTGAT", "AAACATCG"]
    >>> mismatch_dict = get_mismatch_dict(seq_list)
    >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
    True
    """
    mismatch_dict = {}
    for seq in seq_list:
        seq = seq.strip()
        if seq == '':
            continue
        for mismatch_seq in findall_mismatch(seq, n_mismatch):
            mismatch_dict[mismatch_seq] = seq
    return mismatch_dict


def read_one_col(fn):
    """read one column file into list"""
    with open(fn) as f:
        return [x.strip() for x in f]


def parse_pattern(pattern, allowed="CLUNT"):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [slice(0, 8, None), slice(24, 32, None), slice(48, 56, None)]
    >>> pattern_dict['L']
    [slice(8, 24, None), slice(32, 48, None), slice(56, 57, None)]
    """
    pattern_dict = {}
    p = re.compile(r'([A-Z])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        sys.exit(f'Invalid pattern: {pattern}')
    start = 0
    for x, length in tmp:
        if x not in allowed:
            sys.exit(f'Invalid pattern: {pattern}')
        if x not in pattern_dict:
            pattern_dict[x] = []
        end = start + int(length)
        pattern_dict[x].append(slice(start,end))
        start = end
    return pattern_dict

class Barcode:
    """
    """

    def __init__(self, args, display_title=None):
        pass


    @staticmethod
    def get_seq_str_no_exception(seq, sub_pattern_dict):
        """get subseq with intervals in arr and concatenate"""
        return ''.join([seq[item[0]: item[1]] for item in sub_pattern_dict])

    @staticmethod
    def get_seq_str(seq, sub_pattern_dict):
        """
        Get subseq with intervals in arr and concatenate

        Args:
            seq: str
            sub_pattern_dict: [[0, 8], [24, 32], [48, 56]]

        Returns:
            str
            if sequence length is not enough, return ""

        >>> sub_pattern_dict = [[0, 8]]
        >>> seq = "A" * 7
        >>> Barcode.get_seq_str(seq, sub_pattern_dict)
        ""
        >>> seq = "A" * 8
        >>> Barcode.get_seq_str(seq, sub_pattern_dict)
        'AAAAAAAA'
        """
        seq_len = len(seq)
        ans = []
        for item in sub_pattern_dict:
            start, end = item[0], item[1]
            if end > seq_len:
                raise IndexError(f"sequence length is not enough in R1 read: {seq}")
            else:
                ans.append(seq[start:end])
        return ''.join(ans)

    @staticmethod
    def get_seq_list(seq, pattern_dict, abbr):
        """
        >>> pattern_dict = Barcode.parse_pattern("C2L3C2")
        >>> seq = "AAGGGTT"
        >>> Barcode.get_seq_list(seq, pattern_dict, "C")
        ['AA', 'TT']
        """
        
        return [seq[item[0]: item[1]] for item in pattern_dict[abbr]]

    @staticmethod
    def parse_pattern(pattern):
        """
        >>> pattern_dict = Barcode.parse_pattern("C8L16C8L16C8L1U12T18")
        >>> pattern_dict['C']
        [[0, 8], [24, 32], [48, 56]]
        >>> pattern_dict['L']
        [[8, 24], [32, 48], [56, 57]]
        """
        pattern_dict = defaultdict(list)
        p = re.compile(r'([CLUNT])(\d+)')
        tmp = p.findall(pattern)
        if not tmp:
            sys.exit()(f'Invalid pattern: {pattern}')
        start = 0
        for item in tmp:
            end = start + int(item[1])
            pattern_dict[item[0]].append([start, end])
            start = end
        return pattern_dict

    @staticmethod
    def get_abbr_len(pattern_dict, abbr):
        """
        >>> pattern_dict = Barcode.parse_pattern("C8L16C8L16C8L1U12T18")
        >>> Barcode.get_abbr_len(pattern_dict, 'C')
        24
        >>> Barcode.get_abbr_len(pattern_dict, 'L')
        33
        """
        length = 0
        for item in pattern_dict[abbr]:
            length += item[1] - item[0]

        return length
        
    @staticmethod
    def ord2chr(q, offset=33):
        return chr(int(q) + offset)

    @staticmethod
    def qual_int(char, offset=33):
        return ord(char) - offset

    @staticmethod
    def low_qual(quals, minQ, num):
        # print(ord('/')-33)           14
        return True if len([q for q in quals if Barcode.qual_int(q) < minQ]) > num else False

    @staticmethod
    def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
        """
        choose locations where there's going to be a mismatch using combinations
        and then construct all satisfying lists using product

        Return:
        all mismatch <= n_mismatch set. 

        >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
        >>> seq_set = Barcode.findall_mismatch("ACG")
        >>> seq_set == answer
        True
        """
        seq_set = set()
        seq_len = len(seq)
        if n_mismatch > seq_len:
            n_mismatch = seq_len
        for locs in itertools.combinations(range(seq_len), n_mismatch):
            seq_locs = [[base] for base in seq]
            for loc in locs:
                seq_locs[loc] = list(bases)
            for poss in itertools.product(*seq_locs):
                seq_set.add(''.join(poss))
        return seq_set

    @staticmethod
    def get_mismatch_dict(seq_list, n_mismatch=1):
        """
        Return:
        mismatch dict. Key: mismatch seq, value: seq in seq_list

        >>> seq_list = ["AACGTGAT", "AAACATCG"]
        >>> mismatch_dict = Barcode.get_mismatch_dict(seq_list)
        >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
        True
        """
        mismatch_dict = {}

        for seq in seq_list:
            seq = seq.strip()
            if seq == '':
                continue
            for mismatch_seq in Barcode.findall_mismatch(seq, n_mismatch):
                mismatch_dict[mismatch_seq] = seq

        return mismatch_dict

    @staticmethod
    def check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list):
        '''
        Return bool_valid, bool_corrected, corrected_seq

        >>> seq_list = ['ATA', 'AAT', 'ATA']
        >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
        >>> mismatch_dict_list = [Barcode.get_mismatch_dict(['AAA'])] * 3

        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, True, 'AAA_AAA_AAA')

        >>> seq_list = ['AAA', 'AAA', 'AAA']
        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, False, 'AAA_AAA_AAA')
        '''
        bool_valid = True
        bool_corrected = False
        corrected_seq_list = []
        for index, seq in enumerate(seq_list):
            if seq not in correct_set_list[index]:
                if seq not in mismatch_dict_list[index]:
                    bool_valid = False
                    return bool_valid, bool_corrected, ""
                else:
                    bool_corrected = True
                    corrected_seq_list.append(mismatch_dict_list[index][seq])
            else:
                corrected_seq_list.append(seq)

        return bool_valid, bool_corrected, '_'.join(corrected_seq_list)

    @staticmethod
    def parse_whitelist_file(files: list, n_pattern: int, n_mismatch: int):
        """
        files: file paths
        n_pattern: number of sections in pattern
        n_mismatch: allowed number of mismatch bases
        Returns:
            white_set_list
            mismatch_list
        """
        n_files = len(files)
        if n_files == 1 and n_pattern > 1:
            files = [files[0]] * n_pattern
        elif n_files != n_pattern:
            sys.exit(f'number of whitelist files({n_files} files:{files}) != n_pattern({n_pattern})')
        
        white_set_list, mismatch_list = [], []
        for f in files:
            barcodes, _ = read_one_col(f)
            white_set_list.append(set(barcodes))
            barcode_mismatch_dict = Barcode.get_mismatch_dict(barcodes, n_mismatch)
            mismatch_list.append(barcode_mismatch_dict)

        return white_set_list, mismatch_list