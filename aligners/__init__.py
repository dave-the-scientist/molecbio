""" This package implements different sequence alignment algorithms. Currently
it defines the Needleman-Wunsch as the Needleman class.
"""
try:
    from nw_aligner import Needleman, PyNeedleman
except:
    Needleman, PyNeedleman = None, None
