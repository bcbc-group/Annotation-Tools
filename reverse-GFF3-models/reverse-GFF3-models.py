#!/usr/bin/env python
import argparse
import re
from sys import stdout, stderr
from pathlib import Path
from collections import OrderedDict

def regex_parse(text):
    return re.compile(text.encode().decode('unicode_escape'))

parser = argparse.ArgumentParser(
    description='This script reverses gene model names in GFF3 files based on chromosome. It assumes the entries are sorted by postion. It requires the eigth column of each GFF3 row to contain an `ID=*` entry. It also requires that the contig names end with their number.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", help="Path of input GFF3 file.", type=Path)
parser.add_argument(
    "-o", "--output",
    help="Path of output GFF3 file. Prints to STDOUT if not provided.",
    type=Path, default=None, required=False)
parser.add_argument(
    "-m", "--output-map",
    help="Path of output TSV file of old_id/new_id pairs.",
    type=Path, default=None, required=False)
parser.add_argument(
    "-f", "--force",
    help="Clobber existing files.", action='store_true')
parser.add_argument(
    "-p", "--pattern",
    help="pattern describing what part of the ID to swap",
    type=regex_parse, default=re.compile("^[^\.]*"), required=False)
parser.add_argument(
    "-c", "--contig-omit-pattern",
    help="if a contig name matches this regex, it will not be renamed.",
    type=regex_parse, default=re.compile("^Contig.*"), required=False)
args = parser.parse_args()

assert args.input.exists()
assert args.force or args.output is None or not args.output.exists()
assert args.force or args.output_map is None or not args.output_map.exists()


def main():
    contig_ids = OrderedDict()
    with args.input.open() as orig:
        for line in orig:
            if line.startswith("#"): continue
            if line.startswith("##FASTA"): break
            cols = line.strip().split("\t")
            if len(cols) < 8: continue  # only interested in GFF lines
            contig = cols[0]
            if contig in contig_ids:
                if contig_ids[contig]==False:
                    continue
            else:
                omitmatch = args.contig_omit_pattern.match(contig)
                if omitmatch!=None:
                    contig_ids[contig] = False
                    continue
                contig_ids[contig] = OrderedDict();
            old_id_full = getKV(cols[8], "ID")
            old_id = args.pattern.search(old_id_full).group(0)
            if old_id in contig_ids[contig]:
                continue
            else:
                contig_ids[contig][old_id] = True
                
        id_map = OrderedDict();
        for contig in contig_ids:
            if contig_ids[contig]!=False:
                ids = [key for key in contig_ids[contig]]
                id_map.update({key:val for key,val in zip(ids,reversed(ids))})
        if args.output_map!=None:
            with args.output_map.open("w") as output_map:
                for key in id_map:
                    output_map.write(key+"\t"+id_map[key]+"\n")
        
        orig.seek(0) 
        
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
                contig = cols[0]
                if contig_ids[contig]!=False:
                    KVstrings = cols[8].strip().rstrip(";").split(";")
                    KVpairs = [s.split("=") for s in KVstrings]
                    newID = "."
                    for kv in KVpairs:
                        if kv[0] in ("ID","Parent"):
                            old_id = args.pattern.search(kv[1]).group(0)
                            try:
                                kv[1] = re.sub(args.pattern,id_map[old_id],kv[1])
                                newID = kv[1]
                            except KeyError:
                                print(
                                    "\tFailed to process line (name '%s' not found):"
                                    % oldID,
                                    file=stderr)
                                print("\t\t" + line.strip(), file=stderr)
                                
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
