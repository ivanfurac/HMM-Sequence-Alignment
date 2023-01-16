
from itertools import combinations
import random
from random import sample

"""
how many seqs will be considered (max is 19 700)

note: complexity is over O(n choose 2) so be careful
"""

NO_OF_SEQUENCES = 20

def generate_testing_data(pairs):
    """
    generates testing data
    Keeps the aligned version in "testing_data_aligned" directory
    Generates non-aligned pairs by removing all dashes and puts them into "testing_data" directory
    1 pair = 1 file
    pairs are separated with newline character
    """
    i = 0

    for pair in pairs:
        pair = remove_dash_dash(pair[0], pair[1])
        pair = remove_gapx_gapy(pair[0], pair[1])

        with open("testing_data_aligned/aligned_no_"+str(i), "w") as f:
            f.write(pair[0]+"\n"+pair[1])
        
        #remove ALL dashes
        pair = pair[0].replace("-", ""),  pair[1].replace("-", "")

        with open("testing_data/test_no_"+str(i), "w") as f:
            f.write(pair[0]+"\n"+pair[1])

        i += 1

    return

def remove_dash_dash(first: str, second: str):
    """
    removes dashes that occur in both strings at the same place 
    (bc that wouldn't happen in pairwise alignment)

    eg. 
    AC--T
    AC-GA

    becomes
    AC-T
    ACGA
    """
    first_res = ""
    second_res = ""

    for i in range(len(first)):
        if(first[i] == "-" and second[i] == "-"):
            continue
        else:
            first_res += first[i]
            second_res += second[i]

    return first_res, second_res

def remove_gapx_gapy(first, second):
    """
    removes gapx then gapy (or vice versa) situactions
    (bc that wouldn't happen in pairwise alignment)

    eg.
    CA-G
    C-TA

    becomes

    CAG
    CTA
    """
    first_res = ""
    second_res = ""
    i = 0

    while(i < len(first)):

        if(i == len(first)-1): #ovo je jako ružno umorna sam iskrneo, ugl rješava nam error index out of range
            first_res += first[i]
            second_res += second[i]

        elif(first[i] == "-" and second[i+1] == "-"):
            first_res += first[i+1]
            second_res += second[i]
            i += 1

        elif(first[i+1] == "-" and second[i] == "-"):
            first_res += first[i]
            second_res += second[i+1]
            i += 1
        else:
            first_res += first[i]
            second_res += second[i]


        i += 1
      
    return first_res, second_res

def estimate_pi(pairs):
    """
    estimates pi matrix (probability of starting at a particular hidden state)

    estimation by frequency: get state of each pair and return frequency (times_occured/no_of_pairs)

    eg.
    AA...
    C-... 
    starts with mis_match, (can be both match or mismatch hence mis_match)

    -A...
    C-...
    starts with gapx
    """
    pi = {
        "mis_match" : 0,
        "gapy" : 0,
        "gapx" : 0
    }

    for pair in pairs:
        x, y = pair

        if(x[0] == "-" and y[0] != "-"):
            pi["gapx"] += 1
        elif(x[0] != "-" and y[0] == "-"):
            pi["gapy"] += 1
        elif(x[0] != "-" and y[0] != "-"):
            pi["mis_match"] += 1

    n = len(pairs)

    for state in pi.keys():
        pi[state] = pi[state] / n

    return pi

def estimate_A(pairs):
    """
    estimates A* matrix (probability of transitions for hidden states)

    estimation by frequency: find out current state and the next one and +1 that transition (times_this_exact_transition_occured/times_in_that_state_to_any)

    """
    A = {
        "mis_match" : {
            "mis_match" : 0,
            "gapy" : 0,
            "gapx" : 0,
            "end" : 0
        },
        "gapy" : {
            "mis_match" : 0,
            "gapy" : 0,
            "gapx" : 0,
            "end" : 0
        },
        "gapx" : {
            "mis_match" : 0,
            "gapy" : 0,
            "gapx" : 0,
            "end" : 0
        }
    }

    

    for pair in pairs:
        x, y = pair
        i = 0

        current_state = ""
        if(x[0] != "-" and y[0] != "-"):
            current_state = "mis_match"
        elif(x[0] == "-"):
            current_state = "gapx"
        elif(y[0] == "-"):
            current_state = "gapy"
        else:
            ValueError("Unknown state")



        while (i < len(x)-1):
            next_state = ""
            if(x[i+1] != "-" and y[i+1] != "-"):
                next_state = "mis_match"
            elif(x[i+1] == "-"):
                next_state = "gapx"
            elif(y[i+1] == "-"):
                next_state = "gapy"
            else:
                ValueError("Unknown state")

            A[current_state][next_state] += 1
            current_state = next_state
            i += 1

        A[current_state]["end"] += 1

    
    for state in A.keys():
        
        n = sum(A[state].values())

        for next_state in A[state].keys():
            A[state][next_state] /= n


    return A

def estimate_E(pairs):
    """
    estimates E* matrix (probability of emissions for hidden states)

    find out current state and find out current emission and +1 that occurence
    """
    E = {
        "mis_match": {
            "AA": 0,
            "AC": 0,
            "AG": 0,
            "AT": 0,
            "CA": 0,
            "CC": 0,
            "CG": 0,
            "CT": 0,
            "GA": 0,
            "GC": 0,
            "GG": 0,
            "GT": 0,
            "TA": 0,
            "TC": 0,
            "TG": 0,
            "TT": 0,
        },

        "gapy": {
            "A-": 0,
            "C-": 0,
            "G-": 0,
            "T-": 0,
        },

        "gapx": {
            "-A": 0,
            "-C": 0,
            "-G": 0,
            "-T": 0,
        }
    }

    for pair in pairs:
        x, y = pair

        for i in range(len(x)):
            current_state = ""

            if(x[i] != "-" and y[i] != "-"):
                current_state = "mis_match"
            elif(x[i] == "-"):
                current_state = "gapx"
            elif(y[i] == "-"):
                current_state = "gapy"

            try:
                E[current_state][x[i]+y[i]] += 1
            except:
                #E[current_state][x[i]+y[i]] = 1
                print("Key error occured. State " + x[i]+y[i] + " shouldn't exist. Ignored.") 


    for state in E.keys():
    
        n = sum(E[state].values())

        for emission in E[state].keys():
            E[state][emission] /= n
    
    return E


lines = []

with open("hiv_alignment2.fasta", 'r') as f:
    lines = f.read().split(">")

#prvi se mora maknut jer je prazan
lines = sample(lines[1:],NO_OF_SEQUENCES)

#remove all headers
for i in range(0, len(lines)):
    lines[i] = lines[i][lines[i].find("\n") + 1:]
for i in range(len(lines)):
    lines[i] = lines[i].replace("\n", "")

#lines without funny letters 
#only A,C,G,T,- are accepted
clean_lines = [] 
for line in lines:
    clean_lines.append(line)

    for char in line:
        if(char not in ["A", "C", "G", "T", "-"]):
            clean_lines.pop()
            print(char + " found in sequence. Am ignoring it.")
            break

lines = clean_lines
random.shuffle(lines)

#generiraj sve parove
pairs = list(combinations(lines, 2))
random.shuffle(pairs)

print("No of pairs: " + str(pairs))

# 2/3 u training data, 1/3 u testing data
training_pairs = pairs[:int(len(pairs)*2/3)]
testing_pairs = pairs [int(len(pairs)*2/3):]

#očisti crtice na training data
for i in range (0,len(training_pairs)):
    pair = training_pairs[i]
    pair = remove_dash_dash(pair[0], pair[1])
    pair = remove_gapx_gapy(pair[0], pair[1])
    training_pairs[i] = pair


print("Estimating pi...")
pi = estimate_pi(training_pairs)
print(pi)
print("-------------------------------")

print("Estimating A*...")
a = estimate_A(training_pairs)
print(a)
print("-------------------------------")

print("Estimating E*...")
e = estimate_E(training_pairs)
print(e)
print("-------------------------------")

print("Generating training data...")
i = 0
for pair in training_pairs:
    with open("training_data/training_no_"+str(i), "w") as f:
        f.write(pair[0]+"\n"+pair[1])
    i += 1

print("Done.\nGenerating testing data...")
generate_testing_data(testing_pairs)
print("Done.\nSaving matrices...")

values = ""
keys = ""
for key in pi.keys():
    keys += key + "\n"
    values += str(pi[key]) + "\n"

with open("matrix_pi_values", "w") as f:
    f.write(values[:-1])

with open("matrix_pi_keys", "w") as f:
    f.write(keys[:-1])

values = ""
keys = ""
for key in a.keys():
    for subkey in a.keys():
        keys += key + " " + subkey + "\n"
        values += str(a[key][subkey]) + "\n"

with open("matrix_a_values", "w") as f:
    f.write(values[:-1])

with open("matrix_a_keys", "w") as f:
    f.write(keys[:-1])

values = ""
keys = ""
for key in e.keys():
    for subkey in e[key].keys():
        keys += subkey + "\n"
        values += str(e[key][subkey]) + "\n"

with open("matrix_e_values", "w") as f:
    f.write(values[:-1])

with open("matrix_e_keys", "w") as f:
    f.write(keys[:-1])

print("Done.")









