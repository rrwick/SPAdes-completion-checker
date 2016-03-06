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
    output_expected_vs_actual(graph_segments, paths, args.depth, args.min, args.max)


def get_arguments():
    parser = argparse.ArgumentParser(description='SPAdes Completion Checker')
    parser.add_argument('graph', help='The assembly_graph.fastg file made by SPAdes')
    parser.add_argument('paths', help='A file containing one proposed graph path per line')
    parser.add_argument('kmer', type=int, help='The final (largest) kmer used in making this graph')
    parser.add_argument('-d', '--depth', action='store', type=float, help='Read depth cutoff', default=1.0)
    parser.add_argument('--min', action='store', type=float, help='Minimum acceptable ratio')
    parser.add_argument('--max', action='store', type=float, help='Maximum acceptable ratio')

    return parser.parse_args()


def load_graph(filename):
    graph_segments = []
    graph_file = open(filename, 'r')
    for line in graph_file:
        if line.startswith('>'):
            graph_segments.append(GraphSegment(line))
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

def output_expected_vs_actual(graph_segments, paths, depth_cutoff, minimum, maximum):

    print('Segment\t', end="")
    print('\t'.join([str(x) for x in range(1, len(paths) + 1)]), end="")
    print('\tExpected depth\tActual depth\tExpected/Actual')

    positive_graph_segments = [x for x in graph_segments if x.positive]
    positive_graph_segments = sorted(positive_graph_segments, key=lambda x: x.number)

    for segment in positive_graph_segments:

        actual_depth = segment.depth
        if actual_depth < depth_cutoff:
            continue
        number_without_sign = segment.number
        expected_depth = 0.0
        path_occurrences = []
        for path in paths:
            occurrences = path.numbers_without_sign.count(number_without_sign)
            path_occurrences.append(occurrences)
            expected_depth += occurrences * path.depth
        ratio = expected_depth / actual_depth

        if minimum is not None or maximum is not None:
            show = ratio < minimum or ratio > maximum
        else:
            show = True

        if show:
            print(str(number_without_sign) + '\t', end="")
            print('\t'.join([str(x) for x in path_occurrences]), end="")
            print('\t' + str(expected_depth) + '\t' + str(actual_depth) + '\t' + str(ratio))









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
        i = (depth_count - 1) // 2
        if depth_count % 2:
            self.depth = sorted_depths[i]
        else:
            self.depth = (sorted_depths[i] + sorted_depths[i + 1]) / 2.0





if __name__ == '__main__':
    main()