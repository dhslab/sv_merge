from __future__ import division
from __future__ import print_function
import inspect
import sys, argparse, pysam, string, math, re
import pandas as pd
#import pyranges as pr


def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif

def merge_hit(hit1, hit2, outvcf):
    #merge the two hits; keep positions in hit1 because it is the 'primary'
    start = hit1.START1 #min(hit1.START1, hit2.START1)
    stop = hit1.END2 #max(hit1.END2, hit2.END2)
    chrom = hit1.CHROM1
    id_val = f"MERGED:{hit1.INFO1['SVTYPE']}:{hit1.CHROM1}:<{start}:{stop}>"
    filters = (list(set(hit1.FILTERS1 + hit1.FILTERS2))) # + hit2.FILTERS1 + hit2.FILTERS2))) #outvcf.header.filters.keys()
    qual = None #Need to come up with a meaningful merge of the two quality scores
    info = {**hit2.INFO1, **hit1.INFO1}
    if 'ORIGIN' not in info.keys():
        info['ORIGIN'] = hit1.ID1
    if hit1.ID1 not in info['ORIGIN'] and 'MERGED:' not in hit1.ID1:
        info['ORIGIN'] += ',' + hit1.ID1
    if hit2.ID1 not in info['ORIGIN'] and 'MERGED:' not in hit2.ID1:
        info['ORIGIN'] += ',' + hit2.ID1
    format_val = {**hit2.FORMAT1, **hit1.FORMAT1}
    merged_dict = {'start': start, 'stop': stop, 'chrom': chrom, 'id_val': id_val, 'filters': filters, 'qual': qual, 'info': info, 'format_val': format_val, 'ref': 'N', 'alt':hit1.ALT}
    return merged_dict

def update_row(row, hit):
    #updates the row value with information from hit; this is used in merging several SV calls from vcf2 with a single SV call from vcf1
    new_row = pd.DataFrame(index=[row.Index])
    new_row['START1'] = [hit['start']]
    new_row['END1'] = [hit['start']]
    new_row['START2'] = [hit['stop']]
    new_row['END2'] = [hit['stop']]
    new_row['CHROM1'] = [hit['chrom']]
    new_row['CHROM2'] = [hit['chrom']]
    new_row['ID1'] = [hit['id_val']]
    new_row['ID2'] = [hit['id_val']]
    new_row['FILTERS1'] = [hit['filters']]
    new_row['FILTERS2'] = [hit['filters']]
    new_row['QUAL1'] = [hit['qual']]
    new_row['QUAL2'] = [hit['qual']]
    new_row['INFO1'] = [hit['info']]
    new_row['INFO2'] = [hit['info']]
    new_row['FORMAT1'] = [hit['format_val']]
    new_row['FORMAT2'] = [hit['format_val']]
    new_row['REVCOMP'] = row.REVCOMP
    new_row['ORIENT'] = row.ORIENT
    new_row['SVTYPE'] = row.SVTYPE
    new_row['REF'] = row.REF
    new_row['ALT'] = row.ALT
    return tuple(new_row.itertuples())[0]

def vcfdict2vcf(vcfdict, outvcf):
    #Note: the alleles input is a placeholder. Need to implement a solution to retrieve the ref/alt values from the original VCFs.
    vcf = outvcf.new_record(contig=vcfdict['chrom'], start=int(vcfdict['start']), stop=int(vcfdict['stop']), id=vcfdict['id_val'], alleles=(vcfdict['ref'], vcfdict['alt']), filter=vcfdict['filters'], qual=vcfdict['qual'], info=vcfdict['info'])
    for k in vcfdict['format_val'].keys():
        vcf.samples[0][k] = vcfdict['format_val'][k]
    
    return vcf

# function to take a pysam vcf object and return a pandas df in bedpe format
def svtobedpedf(vcf,sample=0):

    verbose = False
    
    # read in vcf records to a bedpe-type pandas object.
    hdr = ['CHROM1','START1','END1','ID1','FILTERS1','QUAL1','INFO1','FORMAT1','CHROM2','START2','END2','ID2','FILTERS2','QUAL2','INFO2','FORMAT2','REVCOMP','ORIENT','SVTYPE','REF','ALT']
    vcfdf = pd.DataFrame(columns = hdr)

    for rec in vcf.fetch():

        # Dragen outputs intervals with no alterations, so we will skip these.
        if rec.alts is None:
            continue       
        
        svtype = ''
        if 'SVTYPE' not in rec.info.keys() or rec.info['SVTYPE'] not in ['DEL','DUP','BND','INV','INS','DUP:TANDEM','CNV']:
            print("***Cannot determine SV type for record: " + str(rec).rstrip(), file=sys.stderr)
            continue

        elif rec.info['SVTYPE'] == 'CNV':
            if len(rec.alts) == 1 and rec.alts[0] == '<DEL>':
                svtype = 'DEL'
                rec.info['SVTYPE'] = 'DEL'
            elif len(rec.alts) == 1 and (rec.alts[0] == '<DUP>' or rec.alts[0] == '<DUP:TANDEM>'):
                svtype = 'DUP'
                rec.info['SVTYPE'] = 'DUP'
            elif (len(rec.alts) == 2 and ('<DEL>' in rec.alts and '<DUP>' in rec.alts)) or 'LOH' in rec.id:
                svtype = 'LOH'
                rec.info['SVTYPE'] = 'LOH'
        else:
            svtype =  rec.info['SVTYPE']

        # get coordinates in bedpe-like format, including slop if there is any
        contig1 = rec.contig
        start1 = rec.pos
        end1 = rec.pos
        if 'CIPOS' in rec.info.keys():
            start1 = max(0, start1 + rec.info['CIPOS'][0]) #the min() call was added to guarantee start1,start2 are non-negative
            end1 = max(0, end1 + rec.info['CIPOS'][1])

        contig2 = rec.contig
        start2 = rec.stop
        end2 = rec.stop
        if rec.stop is not None and 'CIEND' in rec.info.keys():
            start2 = max(0, start2 + rec.info['CIEND'][0]) #the min() call was added to guarantee start1,start2 are non-negative
            end2 = max(0, end2 + rec.info['CIEND'][1])

        # are the SV ends reversed complemented? for dels and dups, no. For inversions, yes. see below for BNDs
        revcomp = 0

        # orientation of strands for BND only. all else None
        orient = None

        if svtype == 'INV':
            revcomp = 1
        
            # for BNDs, need to find the partner locus and the orientation
        if svtype == 'BND':
            # only consider on breakend
            # REF ALT Meaning
            # s t[p[ piece extending to the right of p is joined after t
            # s t]p] reverse comp piece extending left of p is joined after t
            # s ]p]t piece extending to the left of p is joined before t
            # s [p[t reverse comp piece extending right of p is joined before t
            
            m = re.match('([ACGTN]*)(.)(\S+):(\d+)(.)([ACGTN]*)',rec.alts[0])
            contig2 = m.group(3)
            start2 = m.group(4)
            end2 = m.group(4)
        
            if (m.group(1) != '' and m.group(2) == ']') or (m.group(2) == '[' and m.group(6) != ''):
                revcomp = 1

            if (m.group(1) != '' and m.group(2) == '['): # t[p[
                orient = 1
            elif (m.group(1) != '' and m.group(2) == ']'): # t]p]
                orient = 2
            elif (m.group(1) == '' and m.group(2) == ']'): # ]p]t
                orient = 3
            else: # [p[t
                orient = 4

        # now add to df
        info = {}
        for i in rec.info.keys():
            info[i] = rec.info[i]
            
        fmt = {}
        for i in rec.format.keys():
            fmt[i] = rec.samples[sample][i]

        fmt['GT'] = rec.samples[sample].alleles
            
        # if its a BND and the mate has been added, then just update the record, otherwise add it.
        if svtype == 'BND' and 'MATEID' in rec.info.keys() and rec.info['MATEID'][0] in vcfdf['ID1']:
            idx = vcfdf.index[vcfdf['ID1']==rec.info['MATEID']].tolist()
            
            if len(idx) == 1:
                print("too many rows!",file=sys.stderr)
                
            
        else:
            # remember: hdr = ['CHROM1','START1','END1','ID1','FILTERS1','QUAL1','INFO1','FORMAT1','CHROM2','START2','END2','ID2','FILTERS2','QUAL2','INFO2','FORMAT2','REVCOMP','SVTYPE','REF','ALT']            
            lst = [contig1,start1,end1,rec.id,rec.filter.keys(),rec.qual,info,fmt,contig2,start2,end2,rec.id,rec.filter.keys(),rec.qual,info,fmt,revcomp,orient,svtype,rec.ref,','.join(rec.alts)]
            vcfdf = pd.concat([vcfdf, pd.DataFrame.from_records([dict(zip(vcfdf.columns.tolist(),lst))])])

        if verbose:
            print(vcfdf.shape[0],file=sys.stderr)
            
    return vcfdf

#
# Main
#

parser = argparse.ArgumentParser(description='Merge symbolic SV records in two VCFs by SV type with allowed slop/overlap. Records with non-symbolic ALT alleles must match exactly.')
parser.add_argument('vcf1',help='VCF 1, which is the primary VCF')
parser.add_argument('vcf2',help='VCF 2, which is the secondary VCF')
parser.add_argument('outvcf',help='Merged VCF')
parser.add_argument('-n',"--name",help='New sample name to use (default: use primary VCF sample name')
parser.add_argument('-s',"--slop",default=0,help='window around BND breakpoints permitted for merging (default: no slop)')
parser.add_argument('-r',"--roverlap",default=1.0,help='Reciprocal overlap required to merge DEL/DUP (default: 100 percent overlap)')
parser.add_argument('-a',"--all",default=False,help='Include all SV calls in output (default: only merged calls are included)')
parser.add_argument('-i',"--info",default=None,help='Comma-separated INFO tags in secondary VCF to include in merged records. eg: +ABC/ABC will be added, -ABC will be skipped (default: all new tags)')
parser.add_argument('-f',"--format",default=None,help='Comma-separated FORMAT tags in secondary VCF to include in merged records. eg: +ABC/ABC will be added, -ABC will be skipped (default: all new tags; gt field is not added)')
parser.add_argument('-t',"--filters",default=None,help='Filters in secondary VCF to include in merged records. eg: +ABC will be added, -ABC will be skipped (default: all new filters)')
parser.add_argument('-c',"--correct",default=False,help='Adjust symbolic start and end to be the midpoint of the merged records and add IMPRECISE and CIPOS info tags (default: keep primary record coordinates)')

args = parser.parse_args()

sample = str(args.name)
slop = int(args.slop)
roverlap = float(args.roverlap)

if slop and roverlap:
    print("both slop and roverlap entered, roverlap will be used for CNV calls",file=sys.stderr)

cnvtypes = ['DEL','DUP','TANDEM:DUP']
exacttypes = ['INS']
bndtypes = ['BND','INV']

# step 1: load in first VCF into a dict of pyranges objects, keyed by SV type
# step 2: load in second VCF into a dict of pyranges objects, keyed by SV type
# step 3: iterate through first VCF records by SV type. Do overlap with second VCF records and either merge them, or add to output.

vcf1 = pysam.VariantFile(args.vcf1)
vcf2 = pysam.VariantFile(args.vcf2)

# process new tags in vcf2 that should be added/removed
infotags = vcf2.header.info.keys()
if args.info != None:
    tags = str(args.info).split(',')

    # if negated tags are specified, exclude those
    if '-' in str(args.info):
        infotags = Diff(vcf2.header.info.keys(),[it[1:] for it in tags if '-' in it])

    # if any positive tags are specified, only get those
    else:
        infotags = [t for t in vcf2.header.info.keys() if t in [it[1:] for it in tags if '+' in it]]

fmttags = vcf2.header.formats.keys()
if args.format != None:
    tags = str(args.format).split(',')

    # if negated tags are specified, exclude those
    if '-' in str(args.format):
        fmttags = Diff(vcf2.header.formats.keys(),[it[1:] for it in tags if '-' in it])

    # if any positive tags are specified, only get those
    else:
        fmttags = [t for t in vcf2.header.formats.keys() if t in [it[1:] for it in tags if '+' in it]]

filters = vcf2.header.filters.keys()
if args.filters != None:
    tags = str(args.filters).split(',')

    # if negated tags are specified, exclude those
    if '-' in str(args.filters):
        filters = Diff(vcf2.header.filters.keys(),[it[1:] for it in tags if '-' in it])

    # if any positive tags are specified, only get those
    else:
        filters = [t for t in vcf2.header.filters.keys() if t in [it[1:] for it in tags if '+' in it]]
        
# first make output VCF file, copy header, and add new tags and filters
outvcf = pysam.VariantFile(args.outvcf,'w',header=vcf1.header)

# iterate through vcf2 entries and remove records based on user input
for t in vcf2.header.info.keys():
    if t not in infotags:
        vcf2.header.info[t].remove_header()
        
for t in vcf2.header.formats.keys():
    if t not in fmttags:
        vcf2.header.formats[t].remove_header()

for t in vcf2.header.filters.keys():
    if t not in filters:
        vcf2.header.filters[t].remove_header()

# merge the headers
outvcf.header.merge(vcf2.header)

# Adding origin fields to header
info_tag = '##INFO=<ID=ORIGIN,Number=1,Type=String,Description="Origin VCF file.">'
outvcf.header.add_line(info_tag)

vcf1df = svtobedpedf(vcf1)
vcf1df = vcf1df.reset_index(drop=True)
vcf1df = vcf1df.astype({'START1': int, 'START2': int, 'END1': int, 'END2': int})

print("Loaded " + str(vcf1df.shape[0]) + " valid records from " + str(args.vcf1), file=sys.stderr)
vcf1.close()

vcf2df = svtobedpedf(vcf2)
vcf2df = vcf2df.reset_index(drop=True)
vcf2df = vcf2df.astype({'START1': int, 'START2': int, 'END1': int, 'END2': int})

print("Loaded " + str(vcf2df.shape[0]) + " valid records from " + str(args.vcf2), file=sys.stderr)
vcf2.close()

# remember: hdr = ['CHROM1','START1','END1','ID1','FILTERS1','QUAL1','INFO1','FORMAT1','CHROM2','START2','END2','ID2','FILTERS2','QUAL2','INFO2','FORMAT2','REVCOMP','ORIENT','SVTYPE']
count=0

for row in vcf1df.itertuples():
    if row.SVTYPE in cnvtypes: # DEL,DUP,TANDEM:DUP,
        hit = row
        h_ids = [] #these are for keeping track of the indices to be dropped from the df after the merging is complete.
        no_merged = True #Checks whether any variants in vcf2 merged with row
    
        hits = vcf2df[(vcf2df['CHROM1']==row.CHROM1) & \
                      (vcf2df['CHROM2']==row.CHROM2) & \
                      (vcf2df['END2']>=row.START1) & \
                      (vcf2df['START1']<=row.END2)]

        for h in hits.itertuples():
            if row.START1 != row.END1:
                row_left = round((row.START1+row.END1)/2)
            else:
                row_left = row.START1

            if row.START2 != row.END2:
                row_right = round((row.START2+row.END2)/2)
            else:
                row_right = row.END2

            if h.START1 != h.END1:
                h_left = round((h.START1+h.END1)/2)
            else:
                h_left = h.START1

            if h.START2 != h.END2:
                h_right = round((h.START2+h.END2)/2)
            else:
                h_right = h.END2

            overlap = abs(max(row_left, h_left)-min(row_right, h_right))

            if row_left != row_right:
                row_fraction = overlap/abs(row_left-row_right)
            else:
                row_left -= slop
                row_right += slop
                overlap = abs(max(row_left, h_left)-min(row_right, h_right))
                row_fraction = overlap/abs(row_left-row_right)

            if h_left != h_right:
                h_fraction = overlap/abs(h_left-h_right)
            else:
                h_left -= slop
                h_right += slop
                overlap = abs(max(row_left, h_left)-min(row_right, h_right))
                h_fraction = overlap/abs(h_left-h_right)

            if (h_fraction >= roverlap) & (row_fraction >= roverlap):
                hit = merge_hit(hit, h, outvcf)
                row = update_row(row, hit)
                hit = row
                h_ids.append(h.Index)
                no_merged = False

        hit = merge_hit(hit, hit, outvcf)
        if no_merged:
            hit['id_val'] = row.ID1
        hit = vcfdict2vcf(hit, outvcf)
        outvcf.write(hit)
        h_ids = list(set(h_ids))
        vcf1df = vcf1df.drop(row.Index)
        vcf2df = vcf2df.drop(h_ids)

    elif row.SVTYPE in bndtypes: #BND, INV
        hit = row
        h_ids = [] #these are for keeping track of the indices to be dropped from the df after the merging is complete.
        no_merged = True #Checks whether any variants in vcf2 merged with row

        for h in vcf2df[vcf2df['SVTYPE'].isin(bndtypes)].itertuples():
            if row.START1 <= h.START1 and h.START1 - row.END1 <= slop: #slop used for two purposes. Need to fix.
                check1 = True
            elif row.START1 - h.END1 <= slop:
                check1 = True
            if row.START2 <= h.START2 and h.START2 - row.END2 <= slop:
                check2 = True
            elif row.START2 - h.END2 <= slop:
                check2 = True

            if check1 and check2 and row.ORIENT==h.ORIENT:
                hit = merge_hit(hit, h, outvcf)
                row = update_row(row, hit)
                hit = row
                h_ids.append(h.Index)
                no_merged = False

        hit = merge_hit(hit, hit, outvcf)
        if no_merged:
            hit['id_val'] = row.ID1
        hit = vcfdict2vcf(hit, outvcf)
        outvcf.write(hit)
        h_ids = list(set(h_ids))
        vcf1df = vcf1df.drop(row.Index)
        vcf2df = vcf2df.drop(h_ids)

for row in vcf2df.itertuples():
    if row.SVTYPE == 'INS':
        continue
    hit = merge_hit(row, row, outvcf)
    hit['id_val'] = row.ID1
    hit = vcfdict2vcf(hit, outvcf)
    outvcf.write(hit)

for row in vcf1df.itertuples():
    if row.SVTYPE == 'INS':
        continue
    hit = merge_hit(row, row, outvcf)
    hit['id_val'] = row.ID1
    hit = vcfdict2vcf(hit, outvcf)
    outvcf.write(hit)

outvcf.close()
