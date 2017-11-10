
```
usage: rename-GFF3-models.py [-h] [-o OUTPUT] [-f] [-p PREFIX] [-s GENE_START]
                             [-d GENE_COUNT_MULTIPLIER] [-g GENE_DIGITS]
                             [-c CONTIG_DIGITS] [-m MRNA_DIGITS]
                             [-e FEATURE_DIGITS]
                             input

This script renames gene models in GFF3 files. It requires the eigth column of
each GFF3 row to contain an `ID=*` entry, and a `parent=*` entry for any non-
gene annotations. It also requires that the contig names end with their
number.

positional arguments:
  input                 Path of input GFF3 file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path of output GFF3 file. Prints to STDOUT if not
                        provided. 
                        (default: None)
  -f, --force           Clobber existing files. 
                        (default: False)
  -p PREFIX, --prefix PREFIX
                        Prefix for gene IDs. 
                        (default: ANNOT)
  -s GENE_START, --gene-start GENE_START
                        Where to start counting for numbering genes. 
                        (default:5000)
  -d GENE_COUNT_MULTIPLIER, --gene-count-multiplier GENE_COUNT_MULTIPLIER
                        Multiplier for numbering genes. (e.g. gene-count-
                        multiplier=10: 10,20,30) 
                        (default: 10)
  -g GENE_DIGITS, --gene-digits GENE_DIGITS
                        Minimum digits for gene numbering. Numbers with
                        smaller digits are left-padded with zeros. 
                        (default:6)
  -c CONTIG_DIGITS, --contig-digits CONTIG_DIGITS
                        Minimum digits for contig numbering. Numbers with
                        smaller digits are left-padded with zeros. 
                        (default:3)
  -m MRNA_DIGITS, --mrna-digits MRNA_DIGITS
                        Minimum digits for mrna numbering. Numbers with
                        smaller digits are left-padded with zeros. 
                        (default:1)
  -e FEATURE_DIGITS, --feature-digits FEATURE_DIGITS
                        Minimum digits for numbering of other feature types
                        (e.g. exons, cds). Numbers with smaller digits are
                        left-padded with zeros. 
                        (default: 1)
```
