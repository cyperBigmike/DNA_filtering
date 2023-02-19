import filterUtils

datalen = 140
DEPTH = int(input("Which filtering depth between 1-7 is needed? \n"))

if True :
    input_path = "../../Omer-Pilot-Dataset/Design/allseqmerged.fastq"
    output_path = "../../Omer-Pilot-Dataset/output/removeprimers_and_rev_com_primers/reads_filter_without_primers_rev_com_" + str(DEPTH) + ".fastq" 
    output_path_primers = "../../Omer-Pilot-Dataset/output/removeprimers_and_rev_com_primers/reads_filter_with_primers_rev_com_" + str(DEPTH) + ".fastq"
else:
    input_path = "../../Pilot-Nitsan/Design/seq.txt "
    output_path = "../../Pilot-Nitsan/output/removeprimers_and_rev_com_primers/reads_filter_without_primers_rev_com_" + str(DEPTH) + ".fastq" 
    output_path_primers = "../../Pilot-Nitsan/output/removeprimers_and_rev_com_primers/reads_filter_with_primers_rev_com_" + str(DEPTH) + ".fastq"

front_primer = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
back_primer = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

file_input_fastq =open(input_path, "r")
file_output_with_primers = open(output_path_primers, "w")
file_output_without_primers = open(output_path, "w")
lines = file_input_fastq.readlines()


filters = {
    1 : filterUtils.filter_1,
    2 : filterUtils.filter_2,
    3 : filterUtils.filter_3,
    4 : filterUtils.filter_4,
    5 : filterUtils.filter_5,
    6 : filterUtils.filter_6,
    7 : filterUtils.filter_7
}

filters.get(DEPTH, filterUtils.noFilter)(lines, file_output_with_primers, file_output_without_primers,
                                         front_primer, back_primer, datalen)

file_input_fastq.close()
file_output_with_primers.close()
file_output_without_primers.close()

