#!/home/groups/kipac/chto/anaconda2/envs/newenv3/bin/python

from __future__ import division, absolute_import, print_function
import matplotlib
matplotlib.use('Agg')
import os
import sys
import argparse
import redmapper

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calibrate the redMaPPer red sequence')

    parser.add_argument('-c', '--configfile', action='store', type=str, required=True,
                        help='YAML config file')

    args = parser.parse_args()

    calib = redmapper.calibration.RedmapperCalibrator(args.configfile)
    calib.run()

