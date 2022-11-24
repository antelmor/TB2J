#!/usr/bin/env python3
import argparse
import os
import sys
from TB2J.io_merge import merge


def main():
    parser = argparse.ArgumentParser(
        description=
        "TB2J_merge: merge the exchange parameters calculated from different directions. It can accept two types of calculations, first, the lattice structures are rotated from the original structure. The structures could be generated using the TB2J_rotate.py command. The second type consists of different calculations with different spin orientations."
    )
    parser.add_argument(
        'directories',
        type=str,
        nargs='+',
        help=
        'The directories where the TB2J calculations are done. Inside each of them, there should be a TB2J_results directory which contains the magnetic interaction parameter files. e.g. Fe_x Fe_y Fe_z. Alternatively, it could  be the TB2J results directories.'
    )
    parser.add_argument(
        'main_directory',
        type=str,
        default=None,
        help=
        'The directory containning the TB2J calculation of the reference structure. If not given the last element of the "directories" argumetn will be used'
    )
    parser.add_argument(
        '--type',
        '-T',
        type=str,
        default='structure',
        help=
        'The type of calculations, either structure of spin, meaning that the three calculations are done by rotating the structure/spin. '
    )
    args = parser.parse_args()
    merge(*(args.directories), main_path=args.main_directory, method=args.type.strip().lower(), )

main()
