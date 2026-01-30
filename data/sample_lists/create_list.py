#!/usr/bin/env python3

# Output XRootD URLs for a samweb definition
# Program will print to stdout. Redirect to a file using `./create_list.py DEFINITION > file.txt`

import sys
import argparse

import samweb_client


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='create_list',
        description='Output a list of XROOT URLs for files in a SAM definition.')
    parser.add_argument('definition', help='Data set definition name')
    parser.add_argument('-e', '--experiment', help='Experiment variable for SAM', default=None)
    args = parser.parse_args()

    try:
        swclient = samweb_client.SAMWebClient(experiment=args.experiment)
    except samweb_client.client.ExperimentNotDefined:
        print('Must pass "-e EXPERIMENT" option or set environment variable EXPERIMENT, e.g., "export EXPERIMENT=sbnd", to run', file=sys.stderr)
        sys.exit(1)

    try:
        swclient.descDefinition(args.definition)
    except samweb_client.exceptions.DefinitionNotFound:
        print(f'Definition {args.definition} not found.', file=sys.stderr)
        sys.exit(1)

    for f in swclient.listFilesAndLocations(defname=args.definition, schema='root', structured_output=False):
        print(f.split()[0])
