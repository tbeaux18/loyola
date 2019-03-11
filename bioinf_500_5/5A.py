#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 03-10-2019

5A.py

To execute, script utilizes argv
example:

python3 5A.py path-to-file.txt

Input must be in following format

k-integer m-dimension \n
x y \n
x y \n
etc..

Example:
3 2
0.0 0.0
5.0 5.0
0.0 5.0
1.0 1.0
2.0 2.0
3.0 3.0
1.0 2.0
"""

import sys

def handle_data_input(file_input):
    """ takes in the file input, parses, and returns the k center point,
        m dimensions, and a list of tuple'd data points which should equal
        the m dimension number.
        Args:
            file_input (str) : a path to the file to be opened, there are strict
                                rules for its formatting
        Returns:
            k_center (int) : k number provided from input
            data_tuple (lst of tuple) : list of tuples that represent the given
                                        data points
    """

    captured_input = []
    data_points = []

    with open(file_input, 'r') as input_file:
        for line in input_file:
            captured_input.append(line.strip())

    for idx, d_pt in enumerate(captured_input):
        if idx == 0:
            k_m_data = d_pt.split(' ')
            k_center = int(k_m_data[0])
        else:
            data_points.append(d_pt)

    # converts the data points to float objects and tuples them
    data_tuple = [tuple(float(num) for num in data_pair.split(' ')) for data_pair in data_points]

    return k_center, data_tuple



def max_euclideandist_point(data_points, centers):
    """ takes in the data points, and points from centers and returns the
        max point based on the euclidean distance calculation
        Args:
            data_points (lst of tuples) : tuple'd float objects
            centers (lst of tuples) : tuple'd float objects
        Returns:
            max_point
    """
    max_value = 0

    for data_pt in data_points:
        all_distances = []

        for center_pt in centers:
            euclid_dist = 0

            for pts in enumerate(data_pt):
                euclid_dist += (data_pt[pts[0]] - center_pt[pts[0]]) ** 2
            all_distances.append(euclid_dist)

        min_distance = min(all_distances)

        if min_distance > max_value:
            max_value = min_distance
            max_point = data_pt

    return max_point



def farthest_first_traversal(tuple_data, center_k):
    """ takes the list of tuple data, and the k center, randomly generates
        a data point and adds it to the center points set to begin the iteration
        Args:
            tuple_data (lst of tuple) : float objects are tupled in a list
            center_k (int) : the limit of centers
        Returns:
            center_points (set) : set of data points that are considered
                                    the center of each data pair
    """

    center_points = [tuple_data[0]]

    while len(center_points) < center_k:
        center_points.append(max_euclideandist_point(tuple_data, center_points))

    return center_points


def main():
    """ runs main script """

    file_input = sys.argv[1]

    k_center, data_tuple = handle_data_input(file_input)

    k_center_points = farthest_first_traversal(data_tuple, k_center)

    with open('5A-output.txt', 'w') as output:
        for points in k_center_points:
            output.write("{} {}\n".format(points[0], points[1]))

if __name__ == '__main__':
    main()
