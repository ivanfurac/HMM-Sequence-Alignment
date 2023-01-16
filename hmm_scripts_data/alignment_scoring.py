
CONTINUE_MATCH = 5.0
START_MATCH = 5.0
MISMATCH = +1.0
START_GAP = -5.0
CONTINUE_GAP = -1.0

NO_OF_FILES = 40

def get_score(first, second):

    last_state = None
    current_state = None
    score = 0.0

    if(len(first) != len(second)):
        raise Exception()
        print(len(first))
        print(len(second))

    for i in range(0, len(first)):

        x = first[i]
        y = second[i]

        #do the same for gapy
        if(x == "-"):
            current_state = "gapx"
            if(last_state == current_state):
                score += CONTINUE_GAP
            else:
                score += START_GAP

        elif(y == "-"):
            current_state = "gapy"
            if(last_state == current_state):
                score += CONTINUE_GAP
            else:
                score += START_GAP

        elif(x != y):
            current_state = "mismatch"
            score += MISMATCH

        elif(x == y):
            current_state = "match"
            if(last_state == current_state):
                score += CONTINUE_MATCH
            else:
                score += START_MATCH

        last_state = current_state

    return score


all_scores = []
for i in range (0, NO_OF_FILES):
    
    scores = []

    with open("needleman/needleman_no_" + str(i)) as f:
        lines = f.readlines()
        aligned1, aligned2 = lines[0], lines[1]
        aligned1 = aligned1.strip()
        aligned2 = aligned2.strip()[:-1]

        if(len(aligned1) != len(aligned2)):
            print("Unequal legths at needleman/needleman_no_" + str(i) + ": " + str(len(aligned1))+  ", " + str(len(aligned2)))

        scores.append(str(get_score(aligned1, aligned2)))

    with open("testing_data_aligned/aligned_no_" + str(i)) as f:

        lines = f.readlines()
        aligned1, aligned2 = lines[0], lines[1]
        aligned1 = aligned1[:-1]

        if(len(aligned1) != len(aligned2)):
            print("Unequal legths at testing_data_aligned/aligned_no_" + str(i) + ": " + str(len(aligned1)) +", " + str(len(aligned2)))
            print(aligned1[-1])
        scores.append(str(get_score(aligned1, aligned2)))

    #tu ćemo dodat i za naše rezultate

    with open("hmm_aligned/test_no_" + str(i) + ".txt") as f:
        lines = f.readlines()
        aligned1, aligned2 = lines[0], lines[1]

        if(len(aligned1) != len(aligned2)):
            print("Unequal legths at testing_data_aligned/aligned_no_" + str(i) + ": " + str(len(aligned1)) +", " + str(len(aligned2)))

        scores.append(str(get_score(aligned1, aligned2)))

    with open("hmm_rotated/test_no_" + str(i) + ".txt") as f:
        lines = f.readlines()
        aligned1, aligned2 = lines[0], lines[1]

        if(len(aligned1) != len(aligned2)):
            print("Unequal legths at testing_data_rotated/aligned_no_" + str(i) + ": " + str(len(aligned1)) +", " + str(len(aligned2)))

        scores.append(str(get_score(aligned1, aligned2)))


    all_scores.append(scores)



needleman_sum = 0
testing_data_sum = 0
hmm_aligned_sum = 0
hmm_rotated_sum = 0


# scores.tsv file
output_file = "needleman\ttesting_data\thmm_aligned\trotated\n"
for line in all_scores:
    for score in line:
        output_file += str(score) + "\t"

    needleman_sum += float(line[0])
    testing_data_sum += float(line[1])
    hmm_aligned_sum += float(line[2])
    hmm_rotated_sum += float(line[3])

    output_file+="\n"


with open("scores.tsv", "w") as f:
    f.write(output_file)


with open("average_scores", "w") as f:
    f.write("needleman: " + str(needleman_sum/NO_OF_FILES) + "\ntesting_data: " + str(testing_data_sum/NO_OF_FILES) + "\nhmm_aligned: " + str(hmm_aligned_sum/NO_OF_FILES) + "\nhmm_rotated: " + str(hmm_rotated_sum/NO_OF_FILES))
