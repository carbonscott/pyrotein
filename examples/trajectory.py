#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyrotein as pr
import givens as gv
from loaddata import load_xlsx, label_TMs
import colorsimple as cs
from scipy import spatial

# Reverse the order of all elements from element i to element k in array r.
two_opt_swap = lambda r,i,k: np.concatenate((r[:i], r[i:k+1][::-1], r[k+1:]))

def two_opt(cities, improvement_threshold = 0.001,
                    starts_from_first     = False,
                    circular_path         = False):
   """2-opt Algorithm adapted from https://en.wikipedia.org/wiki/2-opt
      if improvement_threshold >= 1, run this many passes.
      if improvement_threshold <  1, exit when improvement is less than this
      much.
   """
   dmat = spatial.distance_matrix(cities, cities)
   if circular_path:
      path_distance = lambda r,c: np.sum([dmat[r[p], r[p-1]]
                                             for p in range(len(r))])
         # Calculate the euclidian distance in n-space of the route r
         # traversing cities c, ending at the path start.
   else:
      path_distance = lambda r: np.sum([dmat[r[p+1], r[p]]
                                           for p in range(len(r)-1)])
         # For a non-circular path (one which ends at a location different
         # from where it starts), use this formula.
   N     = len(cities)
   route = np.arange(N) # np.random.randint(N,size=N)
      # Make an array of row numbers corresponding to cities.
   improvement_factor = 1 # Initialize the improvement factor.
   best_distance = path_distance(route)
      # Calculate the distance of the initial path.
   passes = 0
   while improvement_factor > improvement_threshold and\
                              improvement_threshold < 1 or\
                     passes < improvement_threshold:
      # If the route is still improving, keep going!
      distance_to_beat = best_distance
         # Record the distance at the beginning of the loop.
      for swap_first in range(1 if starts_from_first else 0, len(route)-3):
         # From each city except the first and last,
         for swap_last in range(swap_first+1,len(route)-1):
            # to each of the cities following,
            new_route = two_opt_swap(route,swap_first,swap_last)
               # try reversing the order of these cities
            new_distance = path_distance(new_route)
               # and check the total distance with this modification.
            if new_distance < best_distance:
               # If the path distance is an improvement,
               route = new_route # make this the accepted best route
               best_distance = new_distance
                  # and update the distance corresponding to this route.
      improvement_factor = 1 - best_distance/distance_to_beat
         # Calculate how much the route has improved.
      passes += 1
   return route
      # When the route is no longer improving substantially, stop searching
      # and return the route.


def reverse_sign(u, vh, rank, index_from_zero = True):
    # Comply with the convention (1-based index)
    rank_in_data = rank if index_from_zero else rank - 1

    # Reverse sign...
    u[:, rank_in_data]  = - u[:, rank_in_data]
    vh[rank_in_data, :] = -vh[rank_in_data, :]

    return None


# Specify chains to process...
fl_chain = "chains.comp.xlsx"
lines    = load_xlsx(fl_chain)

# Specify the range of atoms from rhodopsin...
nterm = 1
cterm = 322    # It was 348
backbone = ["N", "CA", "C", "O"]
length_backbone = (cterm - nterm + 1) * len(backbone)

# Load upstream data...
dmats = np.load("dmats.npy")
u     = np.load("u.npy")
s     = np.load("s.npy")
vh    = np.load("vh.npy")

# Allow positive value that fits the distance in nature (optinal)...
reverse_sign(u, vh, 1, index_from_zero = False)
reverse_sign(u, vh, 2, index_from_zero = False)
reverse_sign(u, vh, 4, index_from_zero = False)
reverse_sign(u, vh, 6, index_from_zero = False)

# Calculate the coefficients...
c = np.matmul(np.diag(s), vh)

# Standardize u and c and assign units...
u_ave = np.sqrt(u.shape[0])
c = c / u_ave

# Define a series of rotation...
rotations = [
    [3, 2, -25],
    [3, 4,  10],
    [5, 3,   7],
    [2, 4,   5],
    [5, 2, -10],
    [4, 6,  20],
    [8, 3,   8],
    [6, 5, -40],
    [5, 7,  18],
    [7, 6, -25],
    [7, 2,  20],
    [4, 3,   4],
    [3, 7,  28],
    [7, 6,  10],
    [9, 6,  30],
]
disp_index = -1    # 0-based Python convention
if len(rotations): rank1_last, rank2_last = rotations[disp_index][0:0 + 2]
for rank1, rank2, theta in rotations:
    gv.givens_rotation(u, s, c, rank1, rank2, theta, index_from_zero = False)


# [[[ TSP ]]]
c_ary = np.dstack((c[1], c[2]))[0]
c_dict = {}
id_frst = "6ofj_B"
id_last = "6fkd_A"
id_list = [id_frst]    # Original list of cities
for i, line in enumerate(lines):
    _, pdb, chain = line[:3]

    id = f"{pdb}_{chain}"
    c_dict[id] = c_ary[i]
    if id != id_frst and id != id_last: id_list.append(id)
id_list.append(id_last)

cities = [ c_dict[id] for id in id_list ]
route = two_opt(cities, starts_from_first = True)    # Return index in original cities

fl_export = f'trajectory.dat'
with open(fl_export,'w') as fh:
    for i in route: fh.write(f"{id_list[i]}\n")

new_cities_order = np.array([cities[route[i]] for i in range(len(route))])
with open("test.dat", 'w') as fh:
    for x, y in new_cities_order:
        fh.write(f"{x} {y}\n")
