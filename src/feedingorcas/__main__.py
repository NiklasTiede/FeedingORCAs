""" CLI for feedingORCAs """
# Copyright (c) 2020, Niklas Tiede.
# All rights reserved. Distributed under the MIT License.

import argparse
import logging
import sys


__version__ = "0.1.0"
__author__ = 'Niklas Tiede'
__author_email__ = 'niklastiede2@gmail.com'
__doc_url__ = 'https://feedingorcas.readthedocs.io'
__src_url__ = 'https://github.com/NiklasTiede/feedingORCAs'


def main(argv):
    """
    """
    parser = __build_parser()
    args = parser.parse_args(argv[1:])
    if args.version:
        print(f'feedingORCAs {__version__}')
        return 0


def __build_parser():
    """Constructs the parser for the command line arguments.

    :returns
        An ArgumentParser instance for the CLI.
    """

    parser = argparse.ArgumentParser(description='Cheminformatics toolkit.')
    parser.add_argument('-v',
                        '--version',
                        action='store_true',
                        help='show version number and exit')

    return parser


def run_main():
    try:
        sys.exit(main(sys.argv))
    except Exception as e:
        sys.stderr.write('feedingorcas: ' + str(e) + '\n')
        sys.exit(1)


if __name__ == '__main__':
    run_main()
