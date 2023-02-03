############### HELPER FUNCTIONS ###############
input_path = "/home/seraj/DNAproj/Pilot-Nitsan/seq.txt"
output_path = "/home/seraj/DNAproj/output/removeprimers_and_rev_com_primers/reads_filter_without_primers_rev_com.fastq" 
output_path_primers = "/home/seraj/DNAproj/output/removeprimers_and_rev_com_primers/reads_filter_with_primers_rev_com.fastq"
front_primer = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
back_primer = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
datalen = 140
# front_primer ="TAAGAGACAG"
# back_primer = "CTGTCTCTTA"

def copy_reverse(cp):
    return (cp[::-1])

def copy_compliment(cp):
    cp = cp.rstrip()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(cp)
    bases = [complement[base] for base in bases]
    cp=''.join(bases)
    return (cp)

subsitution_dictionaries={
    'A': {'C', 'G', 'T'},
    'C': {'A', 'G', 'T'},
    'G': {'A', 'C', 'T'},
    'T': {'A', 'C', 'G'}
}

def get_insertion_ball(s): 
    ins_ball = []
    for i in range(len(s)):
        for inserted_char in {'A', 'C', 'G', 'T'}:
            list_s = list(s)
            list_s.insert(i, inserted_char)
            ins_ball.append(''.join(list_s))
    return ins_ball

def get_substitution_ball(s):
    sub_ball = []
    for i in range(len(s)):
        curr_char=s[i]
        for replaced_char in subsitution_dictionaries[curr_char]:
            list_s = list(s)
            list_s[i] = replaced_char
            sub_ball.append(''.join(list_s))
    return sub_ball

def get_del_ball(s):
    del_ball = []
    for i in range(len(s)):
        list_s = list(s)
        del list_s[i]
        del_ball.append(''.join(list_s))
    return del_ball

# input: line(read), front_ind and back_ind(are the indecies from previous checks)
def check_ED(line, front_ind, back_ind, front_ball, back_ball):
    if front_ind == -1:
        for edited_primer in front_ball[0]:
            front_ind = line.find(edited_primer)
            if front_ind != -1:
                break

    if front_ind == -1:
        for edited_primer in front_ball[1]:
            front_ind = line.find(edited_primer)
            if front_ind != -1:
                break

    if front_ind == -1:
        for edited_primer in front_ball[2]:
            front_ind = line.find(edited_primer)
            if front_ind != -1:
                break

    if front_ind == -1:
        return [front_ind, back_ind]
    
    
    if back_ind == -1:
        for edited_primer in back_ball[0]:
            back_ind = line.find(edited_primer)
            if back_ind != -1:
                break

    if back_ind == -1:
        for edited_primer in back_ball[1]:
            back_ind = line.find(edited_primer)
            if back_ind != -1:
                break

    if back_ind == -1:
        for edited_primer in back_ball[2]:
            back_ind = line.find(edited_primer)
            if back_ind != -1:
                break    
    return [front_ind, back_ind]

def aproxemationED(string1, string2):
    if len(string1) > len(string2):
        difference = len(string1) - len(string2)
        # TODO: maybe it's better to consider the insertion happens in the front of
        # the front primer not at the end 
        string1[:difference]

    elif len(string2) > len(string1):
        difference = len(string2) - len(string1)
        string2[:difference]

    else:
        difference = 0

    for i in range(len(string1)):
        if string1[i] != string2[i]:
            difference += 1
    return difference


def aproxemationED_adapt(string1, string2, front):
    if len(string1) > len(string2):
        difference = len(string1) - len(string2)
        # TODO: maybe it's better to consider the insertion happens in the front of
        # the front primer not at the end 
        if front:
            string1[difference:]
        else:    
            string1[:difference]

    elif len(string2) > len(string1):
        difference = len(string2) - len(string1)
        if front:
            string2[difference:]
        else:    
            string2[:difference]

    else:
        difference = 0

    for i in range(len(string1)):
        if string1[i] != string2[i]:
            difference += 1

    return difference
############### **************** ###############

############### IMPLEMENTATIONS  ###############
def find_alone():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    lines = file_input_fastq.readlines() # is this ok with big files?
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")

    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


def find_rev():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    lines = file_input_fastq.readlines()
    # counter = 0
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)

        # if couldn't find, search for revers
        if front_ind == -1 or back_ind == -1:
            line_rev = copy_reverse(line)
            front_ind = line_rev.find(front_primer)
            back_ind = line_rev.find(back_primer)
            line = line_rev
            # if front_ind != -1 and back_ind != -1:
                # counter = counter + 1

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")
    # print(counter)
    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


def find_comp():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    for line in file_input_fastq.readlines():
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)

        # if couldn't find, search for complementary
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            front_ind = line_com.find(front_primer)
            back_ind = line_com.find(back_primer)
            line = line_com

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")
    
    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


def find_comp_rev_OMER():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    # counter = 0
    for line in file_input_fastq.readlines():
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)

        # if couldn't find, search for complementary
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            line_rev_com = copy_reverse(line_com)
            front_ind = line_rev_com.find(front_primer)
            back_ind = line_rev_com.find(back_primer)
            line = line_rev_com
            # if front_ind != -1 and back_ind != -1:
                # counter = 1 + counter

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=10: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")
    # print(counter)
    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()

# the reverse and comp checks here don't add any reads
def find_comp_rev():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    for line in file_input_fastq.readlines():
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)
    #_______________________________________________________#

        # if couldn't find, search for complementary
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            tmp_front_ind = line_com.find(front_primer)
            tmp_back_ind = line_com.find(back_primer)
            if tmp_front_ind != -1 and tmp_back_ind != -1:
                front_ind = tmp_front_ind
                back_ind = tmp_back_ind
                line = line_com

        # if couldn't find, search for revers
        if front_ind == -1 or back_ind == -1:
            line_rev = copy_reverse(line)
            tmp_front_ind = line_rev.find(front_primer)
            tmp_back_ind = line_rev.find(back_primer)
            if tmp_front_ind != -1 and tmp_back_ind != -1:
                front_ind = tmp_front_ind
                back_ind = tmp_back_ind
                line = line_rev

                ### does the order matter? ###
    #_______________________________________________________#

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")
    
    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


# first check reverse complement then ED  
def find_ED():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    
    front_del_ball = get_del_ball(front_primer)
    front_ins_ball = get_insertion_ball(front_primer)
    front_sub_ball = get_substitution_ball(front_primer)
    back_del_ball = get_del_ball(back_primer)
    back_ins_ball = get_insertion_ball(back_primer)
    back_sub_ball = get_substitution_ball(back_primer)
    front_ball = [front_del_ball,front_ins_ball, front_sub_ball]
    back_ball = [back_del_ball, back_ins_ball, back_sub_ball]
    
    lines = file_input_fastq.readlines()
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)
    #_______________________________________________________#
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            line_rev_com = copy_reverse(line_com)
            tmp_front_ind = line_rev_com.find(front_primer)
            tmp_back_ind = line_rev_com.find(back_primer)
            if tmp_front_ind != -1 and tmp_back_ind != -1:
                line = line_rev_com
                front_ind = tmp_front_ind
                back_ind = tmp_back_ind
        
        # if couldn't find, search for primers with ED == 1
        indices = check_ED(line, front_ind, back_ind, front_ball, back_ball)
        front_ind = indices[0]
        back_ind = indices[1]
    #_______________________________________________________#

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")

    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


# first check ED then reverse complement 
def find_comp_rev_oneED():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    
    front_del_ball = get_del_ball(front_primer)
    front_ins_ball = get_insertion_ball(front_primer)
    front_sub_ball = get_substitution_ball(front_primer)
    back_del_ball = get_del_ball(back_primer)
    back_ins_ball = get_insertion_ball(back_primer)
    back_sub_ball = get_substitution_ball(back_primer)
    front_ball = [front_del_ball,front_ins_ball, front_sub_ball]
    back_ball = [back_del_ball, back_ins_ball, back_sub_ball]
    # print(len(back_ball[0]) + len(back_ball[1]) + len(back_ball[2]) )
    
    lines = file_input_fastq.readlines()
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)
    #_______________________________________________________#

        # if couldn't find, search for primers with ED == 1
        indices = check_ED(line, front_ind, back_ind, front_ball, back_ball)
        front_ind = indices[0]
        back_ind = indices[1]
        
        # if couldn't find, search for complementary
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            line_rev_com = copy_reverse(line_com)
            front_ind = line_rev_com.find(front_primer)
            back_ind = line_rev_com.find(back_primer)
            line = line_rev_com
            # if front_ind == -1 or back_ind == -1:
            #     indices = check_ED(line, front_ind, back_ind, front_ball, back_ball)
            #     front_ind = indices[0]
            #     back_ind = indices[1]
    #_______________________________________________________#

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")

    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


def find_comp_rev_aproxemationED():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    
    
    lines = file_input_fastq.readlines()
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)
    #_______________________________________________________#
        if front_ind == -1:
            # if couldn't find, search for primers with ~ED <= 5
            for i in range(len(line) - datalen - len(back_primer) - len(front_primer) + 5) :
                if aproxemationED(line[i : i + len(front_primer)], front_primer) < 6:
                    front_ind = i
        if back_ind == -1 and front_ind != -1:
            for i in range(front_ind + 130, front_ind + 150):
                if aproxemationED(line[i : i + len(back_primer)], back_primer) < 6:
                    back_ind = i
        
        # if couldn't find, search for complementary
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            line_rev_com = copy_reverse(line_com)
            front_ind = line_rev_com.find(front_primer)
            back_ind = line_rev_com.find(back_primer)
            line = line_rev_com
            # if front_ind == -1 or back_ind == -1:
            #     indices = check_ED(line, front_ind, back_ind, front_ball, back_ball)
            #     front_ind = indices[0]
            #     back_ind = indices[1]
    #_______________________________________________________#

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")

    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()


def find_comp_rev_aproxemationED_adapt():
    file_input_fastq =open(input_path, "r")
    file_output_with_primers = open(output_path_primers, "w")
    file_output_without_primers = open(output_path, "w")
    
    
    lines = file_input_fastq.readlines()
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)
    #_______________________________________________________#
    
    # TODO: Needs tuning maybe we can reach more accuracy
        if front_ind == -1:
            # if couldn't find, search for primers with ~ED <= 5
            for i in range(len(line) - datalen - len(back_primer) - len(front_primer) + 5) :
                if aproxemationED_adapt(line[i : i + len(front_primer)], front_primer, True) < 6:
                    front_ind = i
        if back_ind == -1 and front_ind != -1:
            for i in range(front_ind + 130, front_ind + 150):
                if aproxemationED_adapt(line[i : i + len(back_primer)], back_primer, False) < 6:
                    back_ind = i
        
        # if couldn't find, search for complementary
        if front_ind == -1 or back_ind == -1:
            line_com = copy_compliment(line)
            line_rev_com = copy_reverse(line_com)
            front_ind = line_rev_com.find(front_primer)
            back_ind = line_rev_com.find(back_primer)
            line = line_rev_com
            # if front_ind == -1 or back_ind == -1:
            #     indices = check_ED(line, front_ind, back_ind, front_ball, back_ball)
            #     front_ind = indices[0]
            #     back_ind = indices[1]
    #_______________________________________________________#

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")

    file_input_fastq.close()
    file_output_with_primers.close()
    file_output_without_primers.close()

############### **************** ###############


###############  MAIN FUNCTION   ###############


find_alone()
# find_rev() 
# find_comp()
# find_comp_rev() 
# find_comp_rev_OMER()
# find_ED()
# find_comp_rev_oneED()
# find_comp_rev_aproxemationED()
find_comp_rev_aproxemationED_adapt() 
