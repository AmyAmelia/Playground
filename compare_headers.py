#!/usr/bin/env python

# small script to generate cmd line argument to ensure Speed2018 data is in same format as the local SPEED db
# uses the list of Speed db headers  (edited) to format order raw data names

# Read in db headers and add to a list
with open('Speed_col_headers_merge.csv', 'rb') as h:
    line = h.readline()
    temp = str(line, 'utf-8').upper().rstrip("\r\n")
    headers = temp.split(',')
    headers = list(filter(None, headers))
    print("Db headers: n=", len(headers),  headers)
output = open('headers.txt', 'w')
for i in range(len(headers)):
    output.write(headers[i])
    output.write("\n")

# Read in the raw data header and add to a separate list
with open('wgs10k_20180228G-A_speed_pass_rare.meh.all_10lines.txt') as f:
     line1 = f.readline()
    # print(type(line1))
     line1 = line1.upper().rstrip("\r\n")
     look_up = line1.split('\t')
     print("SPEED headers: n=",len(look_up), look_up)


missed = []
cmd = "awk -F $'\\t' '{if (NR!=1) {print "
for i in range(len(headers)):
    if len(cmd) > 34:
        cmd += '"\\t"'

    if headers[i] in look_up:
        indices =look_up.index(headers[i])
        # print(indices +1, look_up[indices])
        cmd = cmd + "$" + str(indices+1)

    else:
        missed.append(headers[i])
        cmd += '"NA"'

        #print (headers[i])
cmd += "}}'"
print(cmd)


# print("missing fields n=", len(missed), missed)