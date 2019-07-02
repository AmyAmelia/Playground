#! /home/amy/anaconda2/bin/python2.7
# Generates an unannotated VCF from the Blat output for DBASS5. header and quality information cloned from an existing vcf and is
# required for later annotation via Inginuity.
# See 'USR' comments for checks/ edits
import csv
import re
from Bio.Seq import Seq

# opens blat output file, and generates a dictionary of outputted variants and alignments,
# adds the name to end of each variant and splits to get Ref and Alt

with open('/home/amy/Documents/blat/170613_DBASS5_blat.tsv', 'rb') as blatout:
    # with open('/home/amy/PycharmProjects/Splicing/output/D5validation_blat.tsv', 'rb') as blatout:
    blat = csv.reader(blatout, delimiter='\t', quotechar='@')
    b = enumerate(blat)
    dBlat = {}
    count = 0
    for number, row in b:
        # skips over blat output heading
        if number > 4 and len(row) > 1:
            # append copy of splice site name to row and split to find ref and alt
            txt = ([i for i in row[9].split('|')])
            # print txt
            # Identify SNVs then splits name to identify ref and alt. Appends to row
            if row[9].find('>') != -1:
                txt1 = ''.join([i for i in txt[1] if not i.isdigit()])
                row.append(filter(None, re.split("[+>-]", txt1)))
                # If variant on +ve strand Add 1 base to end location to get true location of SNV var
                if row[8] == '+':
                    row[16] = str(int(row[16]) + 1)
                # If variant on -ve strand use start position (replace end with start for ease latter on)
                else:
                    row[16] = row[15]
                    # print len(row[21])
                    # print row[21][1], row[21][2]
                    # get reference reverse complement and replace in row
                    refrc = str(Seq(row[21][1]).reverse_complement())
                    row[21][1] = refrc
                    # get Alt reverse complement and replace in row

                    altrc = str(Seq(row[21][2]).reverse_complement())
                    row[21][2] = altrc
                    # print row[21]
                    # print
                n = row[9]

            # Identify ins
            elif row[9].find('ins') != -1:
                txt1 = ''.join([i for i in txt[1]]).replace("+", "+.+")
                row.append(filter(None, re.split("[+]", txt1)))
                # If variant on -ve strand use start position (replace end with start for ease latter on)
                if row[8] == '-':
                    row[16] = row[15]
                    row[21][2] = filter(None, re.split("ins", row[21][2]))
                    # print row

                    # print row[21][2][1]
                # If variant on +ve strand use end position
                else:
                    pass
                n = row[9]

            # Identify dups
            elif row[9].find('dup') != -1:
                # txt = ''.join([i for i in row[9] if not i.isdigit()])
                txt1 = ''.join([i for i in txt[1]]).replace("+", "+.+")
                row.append(filter(None, re.split("[+]", txt1)))
                # If variant on +ve strand Add 1 base to end location to get true location of SNV var
                if row[8] == '+':
                    row[16] = str(int(row[16]) + 1)
                # If variant on -ve strand use start position (replace end with start for ease latter on)
                else:
                    row[16] = row[15]
                n = row[9]
                # print n

            # Identify dels
            elif row[9].find('del') != -1:
                txt1 = ''.join([i for i in txt[1] + "+."])
                row.append(filter(None, re.split("[+]", txt1)))
                # If variant on +ve strand Add 1 base to end location to get true location of SNV var
                if row[8] == '+':
                    row[16] = str(int(row[16]) + 1)
                # If variant on -ve strand use start position (replace end with start for ease latter on)
                else:
                    row[16] = row[15]
                    #print row[21]
                n = row[9]

            else:  # Variants not suitable
                #print row
                n = 'BlackList'

            if n in dBlat:
                dBlat[n].append(row)
            else:
                dBlat[n] = [row]

# csv indexes:
# 0 = match
# 9 = Q name
# 10 = Q size
# 13 = Chr
# 16 = T end

# open existing vcf template, addition information set as variable inorder to use inginuity
ofile = open("./output/DBASS5_Dup.vcf", "ab")
# inf = "\t 2774.01	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-10.832;DP=250;DS;Dels=0.00;FS=1.585;HaplotypeScore=1.9467;" \
#       "MLEAC=1;MLEAF=0.500;MQ=59.48;MQ0=0;MQRankSum=-0.142;QD=11.10;ReadPosRankSum=1.139;SB=-1.367e+03	" \
#       "GT:AD:DP:GQ:PL	0/1:135,115:250:99:2804,0,4183;"
inf ="\t.   PASS    ADP=2293;WT=0;HET=1;HOM=0;NC=0  GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    0/1:255:2293:2293:1090:1203:52.46%:0E0:40:40:170:920:188:1015"
print "variants Blacklisted = ", len(dBlat['BlackList'])
# print dBlat['BlackList'][0]

count = 0
for k in dBlat:
    if k != 'BlackList':  # excludes variants assigned to the blacklist
        # print k
        # print dBlat[k]
        Exmatch = False
        for v in dBlat[k]:  # identify exact matches first and chromosome only vars not haplotypes
            if int(v[0]) == int(v[10]) and len(v[13]) <= 5:
                # print ' > MATCHed 100%'
                try:
                    # print v[13] + "\t" + v[16] + '\t'+ '.' + '\t' + v[21][1] + '\t'+ v[21][2] + '\t'+ v[9]
                    ofile.write(
                        "\n" + v[13] + "\t" + v[16] + '\t' + '.' + '\t' + v[21][1] + '\t' + v[21][2] + inf + ':'+ v[9])
                    count += 1
                    Exmatch = True
                except:
                    print v, "FAILED"
        for v in dBlat[k]:  # if no exact matches found, look of the alignments within 3bp
            if not Exmatch:
                if int(v[0]) >= (int(v[10]) - 3) and len(v[13]) <= 5:
                    # print ' > 3bp of input'
                    try:
                        #     print  v[13] + "\t" + v[16] + '\t' + '.' + '\t' + v[21][1] + '\t' + v[21][2] + v[9]
                        ofile.write(
                            "\n" + v[13] + "\t" + v[16] + '\t' + '.' + '\t' + v[21][1] + '\t' + v[21][2] + inf + ':' + v[9])
                        count += 1
                        Exmatch = True
                    except:
                        print v[9], "FAILED"

    for v in dBlat[k]:
        # if no close matches found, look of the alignments within 20bp (at least 1 var = dup20)
        if not Exmatch:
            if int(v[0]) >= (int(v[10]) - 20) and len(v[13]) <= 5:
                # print ' > 20bp of input'
                try:
                    # print v[13] + "\t" + v[16] + '\t' + '.' + '\t' + v[21][1] + '\t' + v[21][2] + v[9]
                    ofile.write(
                        "\n" + v[13] + "\t" + v[16] + '\t' + '.' + '\t' + v[21][1] + '\t' + v[21][2] + inf + ':' + v[9])
                    count += 1
                    Exmatch = True
                except:
                    print v[9], "FAILED"
                    # print ""

ofile.close()
print count, 'variants added to VCF'
print (len(dBlat) - 1), "unique variants identified "

# remove  duplicates in vcf
lines_seen = set()  # holds lines already seen
outfile = open("./output/final/DBASS5_170925.vcf", "w")
for line in open("./output/DBASS5_Dup.vcf", "r"):
    if line not in lines_seen:  # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()
