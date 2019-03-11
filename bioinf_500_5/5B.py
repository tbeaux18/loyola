#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 03-10-2019

5B.py

To execute, script utilizes argv
example:

python3 5B.py path-to-file.txt

"""

import sys
from math import sqrt

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



def centroid_assignment(data_points, center):

    centroid_distance = []

    for datapt in data_points:

        for centerpt in center:

            dist = round(sqrt((datapt[0] - centerpt[0]) ** 2 \
            + (datapt[1] - centerpt[1]) ** 2), 3)

            centroid_distance.append(dist)

    return centroid_distance



def kmeans_clustering(data_points, k_center):

    # initializing center points with the k number of data points
    center_points = [points for idx, points in enumerate(data_points) if idx < k_center]

    distance = centroid_assignment(data_points, center_points)

    print(distance)


def main():
    """ runs main script """

    file_input = sys.argv[1]

    k_center, data_tuple = handle_data_input(file_input)

    kmeans_clustering(data_tuple, k_center)

if __name__ == '__main__':
    main()
