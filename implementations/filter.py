import filterUtils

input_path = "/home/seraj/DNAproj/Pilot-Nitsan/seq.txt"
output_path = "/home/seraj/DNAproj/output/removeprimers_and_rev_com_primers/reads_filter_without_primers_rev_com.fastq" 
output_path_primers = "/home/seraj/DNAproj/output/removeprimers_and_rev_com_primers/reads_filter_with_primers_rev_com.fastq"
front_primer = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
back_primer = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

file_input_fastq =open(input_path, "r")
file_output_with_primers = open(output_path_primers, "w")
file_output_without_primers = open(output_path, "w")

datalen = 140
DEPTH = 1

filters = {
    1 : filterUtils.filter_1,
    2 : filterUtils.filter_2,
    3 : filterUtils.filter_3,
    4 : filterUtils.filter_4,
    5 : filterUtils.filter_5
}

filters.get(DEPTH, filterUtils.noFilter)(file_input_fastq, file_output_with_primers, file_output_without_primers,
                                         front_primer, back_primer, datalen)

file_input_fastq.close()
file_output_with_primers.close()
file_output_without_primers.close()

