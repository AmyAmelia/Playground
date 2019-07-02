'''
A Python script to parse the txt output of extract_suppl_reads.sh (a Samtools view bam script for identifying split reads in a bam).
and returns a csv filtered for samples that have more than a specified number of split reads.
Output file naming and location is derived from the input file.

Usage, run from command line:
python supread2csv.py -i <input file.txt> -r <OPTIONAL: number of reads observed (positive integer). If not supplied Default = 5>
Flags:
-h --help
-i --input
-r --reads

Author: Amy Slater
Date created: 19th June 2019

'''
import sys
import getopt

def dict_of_reads(filein):
    ''' Loop through text file, parsing samples with reads to a dictionary, for filtering. Output passed records to a csv.'''
    # Try to open text file supplied. Print error if file can not be opened
    try:
        file = open (filein, "r")
    except:
        print("Unable to open ", filein, "\nPlease check file path is correct")
        sys.exit(2)
    # initialise variables to capture info
    bam = ''
    read_list = []
    dict_bam = {}
    output = filein.replace('.txt', '.csv')

    # loop through each row of the text file.
    for row in file:
        # skip empty rows
        if len(row) < 2 :
            # print ('empty row')
            pass
        elif row.endswith('.bam\n'):
            # If previous iterations have identified reads, record the sample and read list as a dictionary before continuing
            if bam != "" and len(read_list) > 0:
                dict_bam[bam] = read_list
                # clear read list
                read_list = []
            # Set bam to current sample, remove carriage return at end of line
            sample = row.split('/')
            bam = sample[-1].strip('\n')
            # print (bam)
        else: # row should represent a bam read. Add to a list.
            # test that row is a a bam read. Bam/Sam have 11 tab seperated required columns.
            try:
                fields = row.split('\t')
                if len(fields) >= 11:
                    read_list.append(row)
            except:
                pass
    # Capture last sample in the dict if bam reads identified
    if bam != "" and len(read_list) > 0:
        dict_bam[bam] = read_list

    print ('Output:', output)
    return dict_bam, output


def build_csv(bam_dict, output, reads):
    '''Pass the dictionary of bam and filter for number of reads. Export as csv'''
    # check a dictionary has been generated.
    if len(bam_dict) > 0:
        # generate and write a header to csv file..
        with open(output, "a+") as ofile:
            header= "Sample\tReadName\tFlag\tRname\tRPos\tMAPQ\tCIGAR\tRnext\tPnext\tTlen\tSeq\tQual\n"
            ofile.write(header)
            for key in bam_dict:
                # loop through dict and filter entries for samples that meet or exceed the reads thresholds
                if len(bam_dict[key]) >= reads:
                    for val in bam_dict[key]:
                        string = key + "\t" + val
                        ofile.write(string )
                # if reads do not exceed threshold.
                else:
                    pass
        ofile.close()
    else:
        print("No reads/samples have been identified from the input file ")


def main(argv):
    ''' Main function to run script to pass the input file, and the desired bam reads observed threshold'''
    try:
        opts, args = getopt.getopt(argv,"hi:r:",["input=","reads="])
        # test that a flag has been submitted
        assert (len(opts)> 0), 'Please supply and input file'
        # initialise variables
        input_file = ''
        reads = 5  # Default reads = 5
    except getopt.GetoptError:
        print ('Input file required. Run as: \n'
               'supreads2cvs.py -i <inputfile> -r <int:reads>')
        sys.exit(2)
    for opt, arg in opts:
        # test that a txt file has been assigned to the input_file flag
        if opt == '-h':
            print ('test.py -i <inputfile> -r <int:reads>')
            sys.exit()
        # test that a txt file has been assigned to the input_file flag
        elif opt in ("-i", "--ifile"):
            assert(arg.endswith('.txt')), 'Please ensure a .txt file has been supplied as the input file'
            input_file = arg
        elif opt in ("-r", "--reads"):
            # test reads flag is supplied as an integer
            try:
                reads = int(arg)
                assert (reads >0), 'Reads supplied must be greater than 1, if in doubt do not provide a -r flag. ' \
                                'Script will run with a default of 5 reads'
            except:
                print('Please ensure -r --reads is passed a positive integer only')
                sys.exit(2)
    print ('Parsing file: ', input_file)
    print ('Read filter:', reads)

    # call functions
    get_bams, output =dict_of_reads(input_file)
    build_csv(get_bams, output, reads)


if __name__ == "__main__":
    print("Initialising supreads2csv.py ....")
    main(sys.argv[1:])