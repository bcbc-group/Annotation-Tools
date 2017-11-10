#!/usr/bin/env python
import argparse
import re
from sys import stdout, stderr
from pathlib import Path

parser = argparse.ArgumentParser(
    description='This script renames gene models in GFF3 files. It requires the eigth column of each GFF3 row to contain an `ID=*` entry, and a `parent=*` entry for any non-gene annotations. It also requires that the contig names end with their number.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", help="Path of input GFF3 file.", type=Path)
parser.add_argument(
    "-o", "--output",
    help="Path of output GFF3 file. Prints to STDOUT if not provided.",
    type=Path, default=None, required=False)
parser.add_argument(
    "-f", "--force",
    help="Clobber existing files.", action='store_true')
parser.add_argument(
    "-p", "--prefix",
    help="Prefix for gene IDs.",
    type=str, default="ANNOT", required=False)
parser.add_argument(
    "-s", "--gene-start",
    help="Where to start counting for numbering genes.",
    type=int, default=5000, required=False)
parser.add_argument(
    "-d", "--gene-count-multiplier",
    help="Multiplier for numbering genes. (e.g. gene-count-multiplier=10: 10,20,30)", 
    type=int, default=10, required=False)
parser.add_argument(
    "-g", "--gene-digits",
    help="Minimum digits for gene numbering. Numbers with smaller digits are left-padded with zeros.",
    type=int, default=6, required=False)
parser.add_argument(
    "-c", "--contig-digits",
    help="Minimum digits for contig numbering. Numbers with smaller digits are left-padded with zeros.",
    type=int, default=3, required=False)
parser.add_argument(
    "-m", "--mrna-digits",
    help="Minimum digits for mrna numbering. Numbers with smaller digits are left-padded with zeros.",
    type=int, default=1, required=False)
parser.add_argument(
    "-e", "--feature-digits",
    help="Minimum digits for numbering of other feature types (e.g. exons, cds). Numbers with smaller digits are left-padded with zeros.",
    type=int, default=1, required=False)
args = parser.parse_args()

assert args.input.exists()
assert args.force or args.output is None or not args.output.exists()

id_formats = {
    'gene': "%s{parent:0%sd}g{N:0%sd}" % (args.prefix, args.contig_digits, args.gene_digits),
    'mrna': "{parent}.{N:0%sd}" % args.mrna_digits,
    'else': "{parent}:{type}:{N:0%sd}" % args.feature_digits
}


def main():
    with args.input.open() as orig:
        name_map = {}
        contig_numbers = {}
        renamed = True
        passes = 0
        print("Creating Name Hierarchy", file=stderr)
        while renamed:
            passes += 1
            print("\t Pass %s... " % passes, end="", flush=True, file=stderr)
            renamed = False
            grouped_items = {
            }  #items whith the same parent are grouped and sorted by position before assigning them names
            # loop over file to create a map of oldIDs->newIDs
            for line in orig:
                if line.startswith("#"): continue
                if line.startswith("##FASTA"): break
                cols = line.split("\t")
                if len(cols) < 5: continue  # only interested in GFF lines
                # find the current ID
                oldID = getKV(cols[8], "ID")
                type, start, stop = cols[2].lower(), int(cols[3]), int(cols[4])
                # fixes problem where CDS are all named the same thing.
                # type needs to be checked on retrieval as well.
                if type == "cds":
                    oldID = oldID + "(%s,%s)" % (start, stop)
                if oldID in name_map: continue  # dont redo already named rows
                parent = None
                if type == "gene":
                    parent = cols[0].strip()  # genes are grouped by contig
                    # find and store contig number if not previously determined
                    if parent not in contig_numbers:
                        numberID = re.search("\d+$", parent)
                        if numberID is None:
                            raise TypeError(
                                "Contig \"%s\" does not end with a number" %
                                parent)
                        contig_numbers[parent] = int(numberID.group(0))
                else:
                    parent = getKV(cols[8], "Parent").split(
                        ","
                    )[0]  # every other type is grouped by its annotated parent
                    if parent not in name_map:
                        continue  # if we haven't renamed the parent, skip on this pass

                # group by parent
                if parent not in grouped_items:
                    grouped_items[parent] = []
                # create a sortable tuple which contains all info needed for renaming
                sort_tuple = (start, stop, oldID, type)
                grouped_items[parent].append(sort_tuple)
            orig.seek(0)  # go back to the start of the file

            # sort grouped entries
            for parent in grouped_items:
                counts = {}
                grouped_items[parent].sort()
                for start, stop, oldID, type in grouped_items[parent]:
                    if type not in counts:
                        counts[type] = 0
                    else:
                        counts[type] += 1

                    newID = None
                    if type == "gene":
                        cnum = contig_numbers[parent]
                        gnum = counts[type] * args.gene_count_multiplier + args.gene_start
                        newID = id_formats["gene"].format(parent=cnum, N=gnum)
                    elif type == "mrna":
                        newparent = name_map[parent]
                        newID = id_formats["mrna"].format(
                            parent=newparent, N=counts[type] + 1)
                    else:
                        newparent = name_map[parent].rsplit('.', 1)[0]
                        newID = id_formats["else"].format(
                            type=type, parent=newparent, N=counts[type] + 1)
                    name_map[oldID] = newID
                    if not renamed: renamed = True
            print("Done.", file=stderr)

        # name map created!
        print("Writing Renamed File", file=stderr)
        with args.output.open(
                "w") if not args.output is None else stdout as output:
            output.write("##gff-version 3\n")
            for line in orig:
                if line.startswith("#"): continue
                if line.startswith("##FASTA"): break
                cols = line.strip().split("\t")
                if len(cols) < 8: continue  # only interested in GFF lines
                type, start, stop = cols[2].lower(), int(cols[3]), int(cols[4])
                KVstrings = cols[8].strip().rstrip(";").split(";")
                KVpairs = [s.split("=") for s in KVstrings]
                newID = ""
                for kv in KVpairs:
                    if kv[0] == ("ID"):
                        oldID = kv[1].strip()
                        if type == "cds":
                            oldID = oldID + "(%s,%s)" % (start, stop)
                        try:
                            kv[1] = name_map[oldID]
                            newID = name_map[oldID]
                        except KeyError:
                            print(
                                "\tFailed to process line (name '%s' not found):"
                                % oldID,
                                file=stderr)
                            print("\t\t" + line.strip(), file=stderr)
                            continue
                    if kv[0] == ("Parent"):
                        oldIDs = [s.strip() for s in kv[1].split(",")]
                        for oldID in oldIDs:
                            try:
                                kv[1] = name_map[oldID]
                            except KeyError:
                                print(
                                    "\tFailed to process line (name '%s' not found):"
                                    % oldID,
                                    file=stderr)
                                print("\t\t" + line.strip(), file=stderr)
                                continue
                for kv in KVpairs:
                    if kv[0] == ("Name"):
                        kv[1] = newID
                newKVstring = ";".join("=".join(p) for p in KVpairs) + ";"
                cols[8] = newKVstring
                newLine = "\t".join(cols) + "\n"
                output.write(newLine)


def getKV(infoString, key):
    IDkv = next(kv for kv in infoString.split(";") if kv.startswith(key + "="))
    return IDkv.split("=")[1].strip()


if __name__ == '__main__':
    main()
