#!/usr/bin/env python
"""
SPAdes Completion Checker

Author: Ryan Wick
email: rrwick@gmail.com
"""

from __future__ import print_function
from __future__ import division
import argparse
import collections


def main():
    args = get_arguments()

    graph_segments = load_graph(args.graph)
    segments_dict = dict([(x.get_number_with_sign(), x) for x in graph_segments])

    paths = load_paths(args.paths, segments_dict)
    determine_graph_segment_uniqueness(graph_segments, segments_dict, paths)
    determine_path_depths(paths, segments_dict, args.kmer)
    output_path_depths(paths, args.out)
    output_table_and_bandage_labels(graph_segments, paths, args.depth, args.min, args.max, args.out)


def get_arguments():
    parser = argparse.ArgumentParser(description='SPAdes Completion Checker')
    parser.add_argument('graph', help='The assembly_graph.fastg file made by SPAdes')
    parser.add_argument('paths', help='A file containing one proposed graph path per line')
    parser.add_argument('kmer', type=int, help='The final (largest) kmer used in making this graph')
    parser.add_argument('min', type=float, help='Minimum acceptable depth ratio')
    parser.add_argument('max', type=float, help='Maximum acceptable depth ratio')
    parser.add_argument('out', help='Prefix for output files')
    parser.add_argument('--depth', action='store', type=float, help='Read depth cutoff', default=1.0)

    return parser.parse_args()


def load_graph(filename):
    graph_segments = []
    graph_file = open(filename, 'r')
    for line in graph_file:
        if line.startswith('>'):
            graph_segments.append(GraphSegment(line))

    # Make a segment 0 to use for path segments that aren't in the graph.
    graph_segments.append(GraphSegment('EDGE_0_length_0_cov_0'))

    return graph_segments


def load_paths(filename, segments_dict):
    paths = []
    paths_file = open(filename, 'r')
    for line in paths_file:
        if line.strip():
            paths.append(ProposedPath(line, segments_dict))
    return paths


def determine_graph_segment_uniqueness(graph_segments, segments_dict, paths):
    for i, path in enumerate(paths):
        other_paths = paths[:i] + paths[i+1:]
        path.determine_unique_segments(other_paths, segments_dict)


def determine_path_depths(paths, segments_dict, kmer):
    for path in paths:
        path.determine_depth(segments_dict, kmer)


def output_path_depths(paths, output_prefix):
    out_file = open(output_prefix + '_path_depths.txt', 'w')
    out_file.write('Path number\tMedian depth\n')
    for i, path in enumerate(paths):
        out_file.write(str(i+1) + '\t' + str(path.depth) + '\n')



def output_table_and_bandage_labels(graph_segments, paths, depth_cutoff, minimum, maximum, output_prefix):

    table_file = open(output_prefix + '_table.txt', 'w')
    bandage_labels_file = open(output_prefix + '_bandage_labels.csv', 'w')

    table_file.write('Segment\t')
    table_file.write('\t'.join(['Path ' + str(x) + ' actual copies\tPath ' + str(x) + ' expected copies' for x in range(1, len(paths) + 1)]))
    table_file.write('\tActual depth\tExpected depth\tActual depth / Expected depth\n')

    bandage_labels_file.write('Segment\t')
    bandage_labels_file.write('\t'.join(['Path ' + str(x) + ' copies: actual, expected' for x in range(1, len(paths) + 1)]))
    bandage_labels_file.write('\tDepth: actual, expected\tDepth: actual/expected\tColour\n')

    positive_graph_segments = [x for x in graph_segments if x.positive]
    positive_graph_segments = sorted(positive_graph_segments, key=lambda x: x.number)

    shown_segments = []

    for segment in positive_graph_segments:

        actual_depth = segment.depth
        if actual_depth < depth_cutoff:
            continue

        # Determine the segment's expected depth, based on the number of occurrences in each path.
        expected_depth = 0.0
        path_occurrences = []
        for path in paths:
            occurrences = path.numbers_without_sign.count(segment.number)
            path_occurrences.append(occurrences)
            expected_depth += occurrences * path.depth
        if expected_depth > 0.0:
            ratio = actual_depth / expected_depth
        else:
            ratio = '-'

        # If the segment only occurs in one path, determine the expected copy number for that path.
        only_one_path = path_occurrences.count(0) == len(paths) - 1
        occurrence_match = False
        if only_one_path:
            expected_occurences = []
            for i, path in enumerate(paths):
                if segment.number in path.numbers_without_sign:
                    expected_occurence = segment.depth / path.depth
                    expected_occurences.append(expected_occurence)
                    occurrence_match = int(round(expected_occurence)) == path_occurrences[i]
                else:
                    expected_occurences.append('-')
        else:
            expected_occurences = ['-'] * len(paths)

        # Decide whether or not to output to the table of results.
        output = True
        if only_one_path and occurrence_match:
            output = False
        if ratio != '-' and minimum is not None or maximum is not None:
            if ratio > minimum and ratio < maximum:
                output = False

        # Write the node to the results table, if appropriate.
        if output:
            table_file.write(str(segment.number))
            for i, occurrences in enumerate(path_occurrences):
                table_file.write('\t' + str(occurrences) + '\t' + str(expected_occurences[i]))
            table_file.write('\t' + str(actual_depth) + '\t' + str(expected_depth) + '\t' + str(ratio) + '\n')
            shown_segments.append(str(segment.number))

        # All nodes are written into the Bandage labels file, regardless of whether they are in the results table.
        bandage_labels_file.write(str(segment.number))
        for i, occurrences in enumerate(path_occurrences):
            rounded_expected_occurence = expected_occurences[i]
            if expected_occurences[i] != '-':
                rounded_expected_occurence = "{0:.2f}".format(rounded_expected_occurence)
            bandage_labels_file.write('\t' + str(occurrences) + ', ' + rounded_expected_occurence)
        rounded_actual_depth = "{0:.1f}".format(actual_depth)
        rounded_expected_depth = "{0:.1f}".format(expected_depth)
        rounded_ratio = ratio
        if ratio != '-':
            rounded_ratio = "{0:.2f}".format(rounded_ratio)
        bandage_labels_file.write('\t' + rounded_actual_depth + ', ' + rounded_expected_depth)
        bandage_labels_file.write('\t' + rounded_ratio + '\t' + get_colour(ratio) + '\n')

    print(', '.join(shown_segments))


def get_colour(ratio):
    if ratio == '-':
        red = 0
        green = 0
        blue = 0
    elif ratio == 1.0:
        red = 127
        green = 127
        blue = 127
    elif ratio < 1.0:
        blueness = 1.0 - ratio
        red = int(127 - blueness * 127)
        green = int(127 - blueness * 127)
        blue = int(127 + blueness * 128)
    elif ratio > 1.0:
        redness = min((ratio - 1.0) / 2.0, 1.0)
        red = int(127 + redness * 128)
        green = int(127 - redness * 127)
        blue = int(127 - redness * 127)
    red = sorted([0, red, 255])[1]
    green = sorted([0, green, 255])[1]
    blue = sorted([0, blue, 255])[1]
    colour_string = '#' + ("%0.2X" % red) + ("%0.2X" % green) + ("%0.2X" % blue)
    return colour_string


class GraphSegment:
    def __init__(self, header):
        header = header.rstrip()
        if header.startswith('>'):
            header = header[1:]
        if header.endswith(';'):
            header = header[:-1]
        header = header.split(':')[0]
        if header.endswith("'"):
            header = header[:-1]
            self.positive = False
        else:
            self.positive = True
        header_parts = header.split('_')
        if len(header_parts) < 6:
            raise Exception('Graph header parts missing')
        self.number = int(header_parts[1])
        self.length = int(header_parts[3])
        self.depth = float(header_parts[5])
        self.expected_depth = 0.0

    def get_number_with_sign(self):
        if self.positive:
            sign = '+'
        else:
            sign = '-'
        return str(self.number) + sign

    def __repr__(self):
        return self.get_number_with_sign()


class ProposedPath:
    def __init__(self, path_string, segments_dict):
        path_string = ''.join(path_string.split()) # Remove all whitespace
        if len(path_string) == 0:
            raise Exception('Empty path')
        self.numbers_with_sign = path_string.split(',')

        # If any of the segments aren't in the graph, change their number to 0.
        fixed_numbers_with_sign = []
        for number in self.numbers_with_sign:
            if number in segments_dict:
                fixed_numbers_with_sign.append(number)
            else:
                fixed_numbers_with_sign.append('0+')
        self.numbers_with_sign = fixed_numbers_with_sign

        self.numbers_without_sign = [int(y.replace('-', '')) for y in [x.replace('+', '') for x in self.numbers_with_sign]]
        self.segments = [segments_dict[x] for x in self.numbers_with_sign]
        self.uniqueness = []
        self.unique_segments = []
        self.depth = 0.0

    def __repr__(self):
        return ','.join([str(x) for x in self.numbers_without_sign])

    def contains_segment(self, number_without_sign):
        return number_without_sign in self.numbers_without_sign

    def determine_unique_segments(self, other_paths, segments_dict):
        uniqueness = []
        for number_without_sign in self.numbers_without_sign:
            unique = True
            for other_path in other_paths:
                if number_without_sign in other_path.numbers_without_sign:
                    unique = False
                    break
            uniqueness.append(unique)
        self.unique_segments = [x for i, x in enumerate(self.numbers_without_sign) if uniqueness[i]]

    def determine_depth(self, segments_dict, kmer):
        depths = []
        for unique_segment in self.unique_segments:
            segment = segments_dict[str(unique_segment) + '+']
            depth = segment.depth
            length = segment.length - kmer
            for i in range(length):
                depths.append(depth)
        sorted_depths = sorted(depths)
        depth_count = len(sorted_depths)
        if depth_count == 0:
            self.depth = 0.0
        else:
            i = (depth_count - 1) // 2
            if depth_count % 2:
                self.depth = sorted_depths[i]
            else:
                self.depth = (sorted_depths[i] + sorted_depths[i + 1]) / 2.0





if __name__ == '__main__':
    main()