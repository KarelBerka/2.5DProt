#!/usr/bin/env python
# -*- coding: UTF-8 -*-


HELIC_MINIMUM = 3
HELIC_PROBABILITY = 0.8
SHEET_MINIMUM = 1
SHEET_PROBABILITY = 0.8
HELIC_MINIMUM2 = 3
HELIC_PROBABILITY2 = 0.0
SHEET_MINIMUM2 = 1
SHEET_PROBABILITY2 = 0.0


def is_helic(name):
    return (name.find("I") != -1) or (name.find("A") != -1) or \
           (name.find("H") != -1) or (name.find("G") != -1) or \
           (name.find("i") != -1) or (name.find("a") != -1) or \
           (name.find("h") != -1) or (name.find("g") != -1)


def is_sheet(name):
    return (name.find("T") != -1) or (name.find("E") != -1) or \
           (name.find("B") != -1) or (name.find("S") != -1) or \
           (name.find("t") != -1) or (name.find("e") != -1) or \
           (name.find("b") != -1) or (name.find("s") != -1)


def is_nontrivial_sse(name, size):
    return (is_helic(name) and size >= HELIC_MINIMUM) or \
           (is_sheet(name) and size >= SHEET_MINIMUM)


ONLY_DISPLAY = [1, 2]
ONLY_DISPLAY_MIN_SIZE = 3


def sse_should_be_shown(sse, h_type, sizes=[]):
    if sse not in h_type:   # not in ONLY_DISPLAY_BY_HTYPE:
        return 1
    if h_type[sse] not in ONLY_DISPLAY:
        return 1
    if sse not in sizes:    # not in ONLY_DISPLAY_BY_HTYPE:
        return 1
    if sizes[sse] < ONLY_DISPLAY_MIN_SIZE:
        return 1
    return 0


def sse_first_level(sse, h_type, sizes=[]):
    if sse not in h_type:   # not in ONLY_DISPLAY_BY_HTYPE:
        return 1
    if h_type[sse] not in [1]:
        return 1
    if sse not in sizes:    # not in ONLY_DISPLAY_BY_HTYPE:
        return 1
    if sizes[sse] < ONLY_DISPLAY_MIN_SIZE:
        return 1
    return 0


def sse_first_and_second_level(sse, h_type, sizes=[]):
    if sse not in h_type:   # not in ONLY_DISPLAY_BY_HTYPE:
        return 1
    if h_type[sse] not in [1, 2]:
        return 1
    if sse not in sizes:    # not in ONLY_DISPLAY_BY_HTYPE:
        return 1
    if sizes[sse] < ONLY_DISPLAY_MIN_SIZE:
        return 1
    return 0


def determine_types(label, prob, size):
    if is_helic(label) and prob > HELIC_PROBABILITY and \
        size >= HELIC_MINIMUM:
        # if it is helix - G,H,I,a and probabily >0.9 and size >4
        out = 1
    elif is_sheet(label) and prob > SHEET_PROBABILITY and \
        size >= SHEET_PROBABILITY:
        # or if it is sheet - T,E,B,S and probabily >0.9 and size >1
        out = 1
    elif is_helic(label) and prob > HELIC_PROBABILITY2 and \
        size >= HELIC_MINIMUM2:
        # if it is helix - G,H,I,a and probabily >0 and size >2
        out = 2
    elif is_sheet(label) and prob > SHEET_PROBABILITY2 and \
        size >= SHEET_MINIMUM2:
        # or if it is sheet - T,E,B,S and probabily >0 and size >1
        out = 2
    else:
        out = 3
    return out
