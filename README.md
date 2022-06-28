# sv_merge
Python tool to merge structural variant calls in two VCF files.

## Requirements
Written in Python 3. Requires [pandas](https://pandas.pydata.org/) and [pysam](https://pysam.readthedocs.io/en/latest/api.html)

## Usage
```bash
python sv_merge.py [options] input1.vcf input2.vcf outfile.vcf
```

Options:
-n, --namehelp #New sample name to use (default: use primary VCF sample name)
-s, --slop #window around BND breakpoints permitted for merging (default: no slop)
-r, --roverlap #Reciprocal overlap required to merge DEL/DUP (default: 100 percent overlap)
-a, --all #Include all SV calls in output (default: only merged calls are included)
-i, --info #Comma-separated INFO tags in secondary VCF to include in merged records. eg: +ABC/ABC will be added, -ABC will be skipped (default: all new tags)
-f, --format #Comma-separated FORMAT tags in secondary VCF to include in merged records. eg: +ABC/ABC will be added, -ABC will be skipped (default: all new tags; gt field is not added)
-t, --filters #Filters in secondary VCF to include in merged records. eg: +ABC will be added, -ABC will be skipped (default: all new filters)
-c, --correct #Adjust symbolic start and end to be the midpoint of the merged records and add IMPRECISE and CIPOS info tags (default: keep primary record coordinates)
