#!/usr/bin/env python3

import sys
import argparse

import samweb_client


SWCLIENT = samweb_client.SAMWebClient()


def definition_exists(defname: str):
    try:
        SWCLIENT.descDefinition(defname)
        return True
    except samweb_client.exceptions.DefinitionNotFound:
        return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='create_list',
        description='Output a list of XROOT URLs for files in a SAM definition.')
    parser.add_argument('definition')
    args = parser.parse_args()

    if not definition_exists(args.definition):
        print(f'Definition {args.definition} not found.', file=sys.stderr)
        sys.exit(1)

    for f in SWCLIENT.listFilesAndLocations(defname=args.definition, schema='root', structured_output=False):
        print(f.split()[0])
