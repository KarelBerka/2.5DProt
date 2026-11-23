#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Functions to work with vectors

import math

def add(first, second):
    assert len(first) == len(second)
    assert len(first) in (2, 3), "This can handle only 2D and 3D vectors now"
    if len(first) == 2:
        return (first[0]+second[0], first[1]+second[1])
    elif len(first) == 3:
        return (first[0]+second[0], first[1]+second[1], first[2]+second[2])

def rm(first, second):
    assert len(first) == len(second)
    assert len(first) in (2, 3), "This can handle only 2D and 3D vectors now"
    if len(first) == 2:
        return (first[0]-second[0], first[1]-second[1])
    elif len(first) == 3:
        return (first[0]-second[0], first[1]-second[1], first[2]-second[2])

def mult(first, num):
    assert len(first) in (2, 3), "This can handle only 2D and 3D vectors now"
    if len(first) == 2:
        return (first[0]*num, first[1]*num)
    elif len(first) == 3:
        return (first[0]*num, first[1]*num, first[2]*num)

def frac(first, num):
    assert len(first) in (2, 3), "This can handle only 2D and 3D vectors now"
    if len(first) == 2:
        return (float(first[0])/num, float(first[1])/num)
    elif len(first) == 3:
        return (float(first[0])/num, float(first[1])/num, float(first[2])/num)


def center_of_line(line):
    """
    Return center of line provided as start & end coords
    """
    normalized = rm(line[1], line[0])
    halfed = frac(normalized, 2)
    center = add(line[0], halfed)
    return center


def angle_target(start, angle, lenght):
    """
    When we know vectors origin, its angle and lenght, what are the coords
    of end of the vector?
    """
    angle_vector = (math.sin(math.radians(angle)), math.cos(math.radians(angle)))
    angle_vector = mult(angle_vector, lenght)
    return add(start, angle_vector)

def vectors_angle(end):
    """
    When we know coords of a vector (start is [0,0]), what is angle it
    have with x axis

    Test:
        vectors_angle((1,1)) => 45
        vectors_angle((-1,1)) => -45
        vectors_angle((-1,-1)) => -135
        vectors_angle((1,-1)) => 135
    """
    return math.degrees(math.atan2(*end))


def vector_length(vector):
    """
    Return length of a vector
    """
    return math.sqrt(vector[0]**2 + vector[1]**2)

def vectors_intersection(start1, end1, start2, end2):
    """
    Return the intersection of two lines (in 2d), defined by
       start1, end1
       start2, end2
    If the intersection is between (start1 and end1) and
    (start2 and end2), then vectors crosses eachother
    """
    x1 = start1[0]
    x2 = end1[0]
    y1 = start1[1]
    y2 = end1[1]
    x3 = start2[0]
    x4 = end2[0]
    y3 = start2[1]
    y4 = end2[1]
    if ((x2-x1)==0) or (((y2-y1)*(x4-x3)-(x2-x1)*(y4-y3))== 0):
       return [None, None]
    else:
       px = ((y2-y1)*(x4-x3)*x1+(x2-x1)*(x4-x3)*(y3-y1)-(y4-y3)*(x2-x1)*x3)/((y2-y1)*(x4-x3)-(x2-x1)*(y4-y3))
       py = (px*(y2-y1)+y1*(x2-x1)-x1*(y2-y1))/(x2-x1)
       return [px, py] 
