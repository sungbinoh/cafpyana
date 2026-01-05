#!/usr/bin/env python3

import re
import sys
import argparse
import pathlib

import samweb_client


SWCLIENT = samweb_client.SAMWebClient()
XROOT_URL = "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr"
XROOT_RE = re.compile(r'^/pnfs(.*)')


def definition_exists(defname: str):
    try:
        SWCLIENT.descDefinition(defname)
        return True
    except samweb_client.exceptions.DefinitionNotFound:
        return False


def p2x(pnfs_path: str, check=False):
    """Convert a pnfs path to xrootd url."""
    if check:
        # do a check, but this is slow
        p = pathlib.Path(pnfs_path)
        if not p.is_file():
            raise RuntimeError(f'File {pnfs_path} not found when converting to XROOT URL')
        if not XROOT_RE.match(pnfs_path):
            raise RuntimeError(f'File {pnfs_path} does not look like a /pnfs path')

    return XROOT_RE.sub(rf'{XROOT_URL}\1', pnfs_path)


def main(defname: str):
    missing_files = []
    iterator = SWCLIENT.listFiles(defname=defname, stream=True)
    for fname, locs in SWCLIENT.locateFilesIterator(iterator, chunksize=1000):
        if not locs:
            missing_files.append(fname)
            continue

        loc = locs[0]['location']
        if not loc.startswith('enstore:/pnfs'):
            missing_files.append(fname)
            continue

        print(p2x(loc.removeprefix('enstore:') + '/' + fname))

    if missing_files:
        print('Missing files', file=sys.stderr)
        for f in missing_files:
            print(f' {f}', file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='create_list',
        description='Output a list of XROOT URLs for files in a SAM definition.')
    parser.add_argument('definition')
    args = parser.parse_args()

    if not definition_exists(args.definition):
        print(f'Definition {args.definition} not found.', file=sys.stderr)
        sys.exit(1)

    main(args.definition)
