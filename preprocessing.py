
from itertools import combinations

"""
how many seqs will be considered (max is 19 700)

note: complexity is over O(n choose 2) so be careful
"""
NO_OF_SEQUENCES = 100


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
def remove_dash_dash(first: str, second: str):
    first_res = ""
    second_res = ""

    for i in range(len(first)):
        if(first[i] == "-" and second[i] == "-"):
            continue
        else:
            first_res += first[i]
            second_res += second[i]

    return first_res, second_res


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
def remove_gapx_gapy(first, second):
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
def estimate_pi(pairs):
    pi = {
        "mis_match" : 0,
        "gapx" : 0,
        "gapy" : 0
    }

    for pair in pairs:
        x, y = pair
        print(x[0] + ", " + y[0])

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

"""
estimates A* matrix (probability of transitions for hidden states)

estimation by frequency: find out current state and the next one and +1 that transition (times_this_exact_transition_occured/times_in_that_state_to_any)

"""
def estimate_A(pairs):
    A = {
        "mis_match" : {
            "mis_match" : 0,
            "gapx" : 0,
            "gapy" : 0,
            "end" : 0
        },
        "gapx" : {
            "mis_match" : 0,
            "gapx" : 0,
            "gapy" : 0,
            "end" : 0
        },
        "gapy" : {
            "mis_match" : 0,
            "gapx" : 0,
            "gapy" : 0,
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


"""
estimates E* matrix (probability of emissions for hidden states)

find out current state and find out current emission and +1 that occurence
"""
def estimate_E(pairs):

    #ovo je bože sačuvaj... i na kraju sam s try except riješila ona luda slova
    E = {
        "mis_match": {
            "AA": 0,
            "AG": 0,
            "AC": 0,
            "AT": 0,
            "GA": 0,
            "GG": 0,
            "GC": 0,
            "GT": 0,
            "CA": 0,
            "CG": 0,
            "CC": 0,
            "CT": 0,
            "TA": 0,
            "TG": 0,
            "TC": 0,
            "TT": 0,
            "TY": 0,
        },

        "gapx": {
            "-A": 0,
            "-C": 0,
            "-G": 0,
            "-T": 0,
            "-Y": 0,
            "-R": 0,
            "-W": 0,
            "-M": 0,
        },

        "gapy": {
            "A-": 0,
            "C-": 0,
            "G-": 0,
            "T-": 0,
            "Y-": 0,
            "R-": 0,
            "W-": 0,
            "M-": 0
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
                E[current_state][x[i]+y[i]] = 1

    for state in E.keys():
    
        n = sum(E[state].values())

        for emission in E[state].keys():
            E[state][emission] /= n
    
    return E


lines = []
len(lines)

with open("hiv_alignment.fasta", 'r') as f:
    lines = f.read().split(">")

#maknuti (prvi se mora maknut jer je prazan)
lines = lines[1:NO_OF_SEQUENCES+1]

#remove all headers
for i in range(0, len(lines)):
    lines[i] = lines[i][lines[i].find("\n") + 1:]
for i in range(len(lines)):
    lines[i] = lines[i].replace("\n", "")



pairs = list(combinations(lines, 2))
print(len(pairs))

for i in range (0,len(pairs)):
    pair = pairs[i]
    pair = remove_dash_dash(pair[0], pair[1])
    pair = remove_gapx_gapy(pair[0], pair[1])
    pairs[i] = pair

print("estimated PI*: ")
pi = estimate_pi(pairs)
print(pi)
print("-------------------------------")

print("estimated A*: ")
a = estimate_A(pairs)
print(a)
print("-------------------------------")

print("estimated E*: ")
e = estimate_E(pairs)
print(e)
print("-------------------------------")


