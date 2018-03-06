
```
usage: reverse-GFF3-models.py [-h] [-o OUTPUT] [-m OUTPUT_MAP] [-f]
                              [-p PATTERN] [-c CONTIG_OMIT_PATTERN]
                              input

This script reverses gene model names in GFF3 files based on chromosome. It
assumes the entries are sorted by postion. It requires the eigth column of
each GFF3 row to contain an `ID=*` entry. It also requires that the contig
names end with their number.

positional arguments:
  input                 Path of input GFF3 file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path of output GFF3 file. Prints to STDOUT if not
                        provided. (default: None)
  -m OUTPUT_MAP, --output-map OUTPUT_MAP
                        Path of output TSV file of old_id/new_id pairs.
                        (default: None)
  -f, --force           Clobber existing files. (default: False)
  -p PATTERN, --pattern PATTERN
                        pattern describing what part of the ID to swap
                        (default: re.compile('^[^\\.]*'))
  -c CONTIG_OMIT_PATTERN, --contig-omit-pattern CONTIG_OMIT_PATTERN
                        if a contig name matches this regex, it will not be
                        renamed. (default: re.compile('^Contig.*'))
```
