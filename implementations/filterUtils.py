ED_approx_bound = 7

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


# considering the size diff as insertions, and the diff between the remainings as substitution
def aproxemationED(string1, string2):
    if len(string1) > len(string2):
        difference = len(string1) - len(string2)
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

                #######################                           #######################
                #######################    FILTERING FUNCTIONS    ####################### 
                #######################                           #######################

# DEPTH = 1  find_alone
def filter_1(file_input_fastq, file_output_with_primers, file_output_without_primers,
             front_primer, back_primer, datalen):

    lines = file_input_fastq.readlines() # is this ok with big files?
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - datalen) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")


# DEPTH = 2  find_comp_rev_OMER
def filter_2(file_input_fastq, file_output_with_primers, file_output_without_primers,
             front_primer, back_primer, datalen):

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
            
            if abs(len(seq_without_primers) - datalen) <= 5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")



# DEPTH = 3   find_comp_rev_aproxemationED
def filter_3(file_input_fastq, file_output_with_primers, file_output_without_primers,
             front_primer, back_primer, datalen):
    
    lines = file_input_fastq.readlines()
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)
    #_______________________________________________________#
        if front_ind == -1:
            # if couldn't find, search for primers with ~ED <= 5
            for i in range(len(line) - datalen - len(back_primer) - len(front_primer) + 5) :
                if aproxemationED(line[i : i + len(front_primer)], front_primer) < ED_approx_bound:
                    front_ind = i
        if back_ind == -1 and front_ind != -1:
            for i in range(front_ind + 130, front_ind + 150):
                if aproxemationED(line[i : i + len(back_primer)], back_primer) < ED_approx_bound:
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
            
            if abs(len(seq_without_primers) - datalen) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")


# DEPTH = 4   find_comp_rev_oneED
def filter_4(file_input_fastq, file_output_with_primers, file_output_without_primers,
             front_primer, back_primer, datalen):
    
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
    #_______________________________________________________#

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - datalen) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")


# DEPTH = 4  find_alone_onePrimer()
def filter_5(file_input_fastq, file_output_with_primers, file_output_without_primers,
             front_primer, back_primer, datalen):
    print("With this filter you might get more and faster reads, however, the accuracy is bad")
    lines = file_input_fastq.readlines()
    for line in lines:
        front_ind = line.find(front_primer)
        back_ind = line.find(back_primer)

        # if we succeeded to find only one of the primers, calculate the possible position of the other
        if front_ind != -1 and back_ind == -1:
            back_ind = front_ind + len(front_primer) + datalen
        
        elif front_ind == -1 and back_ind != -1:
            front_ind = back_ind - datalen - len(front_primer)

        if front_ind != -1 and back_ind != -1:
            seq_with_primers = line[front_ind:back_ind + len(back_primer)]
            seq_without_primers = line[front_ind + len(front_primer): back_ind]
            
            if abs(len(seq_without_primers) - 140) <=5: # we take data in offset 5 at most ,else throw it.
                file_output_with_primers.write(seq_with_primers + "\n")
                file_output_without_primers.write(seq_without_primers + "\n")




def noFilter(input_path, output_path, output_path_primers, front_primer, back_primer, datalen):
    print("This Depth isn't supported")