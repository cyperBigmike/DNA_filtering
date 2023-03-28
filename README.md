# DNA_filtering
filtering DNA sequences , filters levels : 

 1 : basic find
 
 2 : complement reverse 
 
 3 : Hamming distance approximation 
 
 4 : one edit distance
 
 5 : filter 1 - 4 combined
 
 6 : find sub primer
 
 7 : find One edit distnace with sub primer
 
 8 : find one primer



**Usage - 
you need to provide the input file path and what filter level you want to use like this  : 

' python3.6 filter.py  input.fastq 1 '

**output - 

you will get the output files in the current directory like this :

    reads_filter_without_primers_rev_com_1.fastq
    
    reads_filter_with_primers_rev_com_ 1.fastq
