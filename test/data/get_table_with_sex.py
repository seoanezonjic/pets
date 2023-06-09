import random
random.seed(42)

new_table = []

with open("./100_test_dataset.txt") as f:
    header = True
    for line in f:
        if header:
            header_names = line.strip().split("\t") + ["sex"]
            new_table.append(header_names)
            header = False
            continue
        line = line.strip().split("\t")
        line.append(random.choice(["M", "F"]))
        new_table.append(line)

with open("./100_test_dataset_with_sex.txt", "w") as f:
    for line in new_table:
        f.write("\t".join(line) + "\n")