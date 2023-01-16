import os
from minineedle import needle

no_of_pairs = len(os.listdir('testing_data'))

alignment = needle.NeedlemanWunsch("ACA", "C")
alignment.align()
a = alignment.__str__().split("\n") #ne pitaj, vjeruj
print(alignment)
print(a[1][1:] + "\n" + a[2][1:])


for file_index in range(10, 40):
    with open("testing_data/test_no_"+str(file_index)) as f:
        lines = f.readlines()

    first, second = lines[0], lines[1]

    print(file_index)
    alignment = needle.NeedlemanWunsch(first, second)
    alignment.align()


    with open("needleman/needleman_no_"+str(file_index), "w") as f:
        a = alignment.__str__().split("\n") #ne pitaj, vjeruj
        f.writelines(a[1][1:] + "\n" + a[3][1:])
