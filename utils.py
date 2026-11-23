#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math
import numpy as np
import tabulate

const = 4


def count_angle_in_triangle(a, b, c):
    # count angle near vertex b
    aa = np.array([a[0], a[1], a[2]])
    bb = np.array([b[0], b[1], b[2]])
    cc = np.array([c[0], c[1], c[2]])
    ba = aa - bb
    bc = cc - bb
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


def count_angle_in_triangle_2d(a, b, c):
    # count angle near vertex b - in 2D
    aa = np.array([a[0], a[1], 0])
    bb = np.array([b[0], b[1], 0])
    cc = np.array([c[0], c[1], 0])
    ba = aa - bb
    bc = cc - bb
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def vector_lenght(a):
    return math.dist(a, (0, 0, 0))

# from given 2d layout find the most distant pair and returns it
# if there is an input set of points, then the function choose
# only from the given set
def find_most_distant_pair(layout, set):
    max_dest = 0
    pair = ["",""]
    for sse1 in set:
       for sse2 in set:
          dest = math.dist(layout[sse1], layout[sse2])
          if dest > max_dest:
             pair = [sse1, sse2]
             max_dest = dest
    return pair

# from given 2d layout find the most distant point
def find_most_distant_point(layout, set, start):
    max_dest = -1
    end = ""
    for sse in set:
          dest = math.dist(layout[start], layout[sse])
          if dest > max_dest:
             end = sse
             max_dest = dest
    return end


# returns the distance of line l = a+kb and point p
def point_to_line_dist(a,b,p):
   for i in range(const+1):
      point = [a[0]+i/(const*1.0)*b[0], a[1]+i/(const*1.0)*b[1], a[2]+i/(const*1.0)*b[2]]
      if (i == 0):
         mi = math.dist(a,p)
      else:
         m = math.dist(point,p)
         if (m<mi):
            mi = m
   return mi

# returns the distance of line l = s+kd and point p
def point_to_line_dist2(s,d,p):
   aux = (d[1]**2+d[0]**2)
   if aux == 0:
      return 0
   else:
      return abs((d[1]*p[0]-d[0]*p[1]+(s[1]*d[0]-s[0]*d[1])))/math.sqrt(aux)

# returns 0 if the minimal distance from point p to line from s to e
# is reached inside line segment s-e
def point_to_line_min(s,e,p):
   d = [s[0]-e[0],s[1]-e[1]]
   const = (s[1]*d[0]+p[0]*d[1]-s[0]*d[1]-p[1]*d[0])/(d[1]*d[1]+d[0]*d[0])
   xcoord = p[0]-const*d[1]
   ycoord = p[1]+const*d[0]
   if (((s[0]<xcoord) and (xcoord<e[0])) or \
      ((s[0]>xcoord) and (xcoord>e[0]))) and \
      (((s[1]<ycoord) and (ycoord<e[1])) or \
      ((s[1]>ycoord) and (ycoord>e[1]))):
      return 0
   else:
      return 1

# returns the distance of lines l_1 = a+kb, l_2 = c+ld
# return the optimalisation of the distance (can't detect the problem)
# grid const x const test all distances choose the minimum
def line_to_line_dist(a,b,c,d):
   for i in range(const+1):
      point1 = [a[0]+i/(const*1.0)*b[0], a[1]+i/(const*1.0)*b[1], a[2]+i/(const*1.0)*b[2]]
      for j in range(const+1):
         point2 = [c[0]+j/(const*1.0)*d[0], c[1]+j/(const*1.0)*d[1], c[2]+j/(const*1.0)*d[2]]
         if (i == 0) and (j == 0):
            mi = math.dist(a,c)
         else:
            m = math.dist(point1,point2)
            if (m<mi):
               mi = m
   return mi

def rotate_point(centerPoint,point,angle):
    """
    Rotates a 2D point around another centerPoint. Angle is in radians.
    Rotation is counter-clockwise
    See http://stackoverflow.com/questions/20023209/function-for-rotating-2d-objects
    """
    temp_point = (point[0] - centerPoint[0], point[1] - centerPoint[1])
    temp_point = (temp_point[0] * math.cos(angle) - temp_point[1] * math.sin(angle), temp_point[0] * math.sin(angle) + temp_point[1] * math.cos(angle))
    return (temp_point[0] + centerPoint[0], temp_point[1] + centerPoint[1])


def get_distance_matrix4n(vectors):
#  TODO time consuming if I use version 0/2 in new generate script
    """Get mutual minimal distances of line segments aproximating each helic"""
    return get_distance_matrix4n_2(vectors, vectors)


def get_distance_matrix4n_2(vectors1, vectors2):
#  TODO time consuming if I use version 0/2 in new generate script
    """Get mutual minimal distances of line segments aproximating each helic"""
    distances = {}
    for h1 in vectors1:
        if h1 not in distances:
            distances[h1] = {}
        h1_start = vectors1[h1][0]
        h1_end = vectors1[h1][1]
        h1_diff = [h1_end[0]-h1_start[0], h1_end[1]-h1_start[1], h1_end[2]-h1_start[2]]
        for h2 in vectors2:
            h2_start = vectors2[h2][0]
            h2_end = vectors2[h2][1]
            h2_diff = [h2_end[0]-h2_start[0], h2_end[1]-h2_start[1], h2_end[2]-h2_start[2]]
            distances[h1][h2] = line_to_line_dist(h1_start, h1_diff, h2_start, h2_diff)
    return distances

def get_distance_matrix4(protein, chain_id, domain):
#  TODO time consuming if I use version 0/2 in new generate script
  """Get mutual minimal distances of line segments aproximating each helic"""
  distances = {}
  secstr = protein.get_helices(chain_id, domain)
  for h1 in secstr:
    if h1 not in distances:
      distances[h1] = {}
    h1_line = protein.get_helic_aprox_bestline(chain_id, domain, h1)
    h1_diff = h1_line[1] - h1_line[0]
    h1_start = h1_line[0]
    for h2 in secstr:
      h2_line = protein.get_helic_aprox_bestline(chain_id, domain, h2)
      h2_diff = h2_line[1] - h2_line[0]
      h2_start = h2_line[0]
      distances[h1][h2] = line_to_line_dist(h1_start, h1_diff, h2_start, h2_diff)
  return distances

def get_distance_matrix4_2(protein, chain_id1, domain1, chain_id2, domain2, short = False):
  """Get mutual minimal distances of line segments aproximating each helic"""
  distances = {}
  secstr1 = protein.get_helices(chain_id1, domain1)
  secstr2 = protein.get_helices(chain_id2, domain2)
  if short == True:
       # short version uses only compute distance matrix only for 3 sses in given domain
       # it is used for big proteins with more then 70 domains (class Layout) -
       # in this case time complexity is too big and the gains of long version are minimal
       ss1 = {}
       for i in secstr1.keys():
           ss1[i] = secstr1[i]
           if len(ss1) == 3:
               break
       ss2 = {}
       for i in secstr2.keys():
           ss2[i] = secstr2[i]
           if len(ss2) == 3:
               break
       secstr1 = ss1
       secstr2 = ss2
  for h1 in secstr1:
    if h1 not in distances:
      distances[h1] = {}
    h1_line = protein.get_helic_read_line22(chain_id1, domain1, h1)
    h1_diff = np.array(h1_line[1]) - np.array(h1_line[0])
    h1_start = h1_line[0]
    for h2 in secstr2:
      h2_line = protein.get_helic_read_line22(chain_id2, domain2, h2)
      h2_diff = np.array(h2_line[1]) - np.array(h2_line[0])
      h2_start = h2_line[0]
      dist = line_to_line_dist(h1_start, h1_diff, h2_start, h2_diff)
      distances[h1][h2] = dist
  return distances


def get_start_end_distance_matrix(protein, chain_id, domain):
    """Get mutual minimal distances of line segments aproximating each helic"""
    se_distances = {}
    secstr = protein.get_helices(chain_id, domain)
    for h1 in secstr:
        if h1 not in se_distances:
           se_distances[h1] = {}
        h1_line = protein.get_helic_aprox_bestline(chain_id, h1)
        for h2 in secstr:
            h2_line = protein.get_helic_aprox_bestline(chain_id, h2)
            se_distances[h1][h2] = [math.dist(h1_line[0], h2_line[0]), \
                                    math.dist(h1_line[0], h2_line[1]), \
                                    math.dist(h1_line[1], h2_line[0]), \
                                    math.dist(h1_line[1], h2_line[1])]
    return se_distances


def get_table(data, roundme=2, sortme=True, tablefmt='simple'):
  table = []
  header_up = []
  for row in data.values():
    for col in row:
        if col not in header_up:
            header_up.append(col)
  header_left = list(data.keys())
  if sortme:
    header_up.sort()
    header_left.sort()
  for row in header_left:
    table_row = [row]
    for column in header_up:
      try:
        cell = data[row][column]
      except KeyError:
        cell = None
      if roundme >= 0:
        if isinstance(cell, float):
          table_row.append(round(cell * 10**roundme) / 10**roundme)
        else:
          table_row.append(str(cell))
      else:
        table_row.append(cell)
    table.append(table_row)
  return tabulate.tabulate(table, headers=['']+header_up, tablefmt=tablefmt)
