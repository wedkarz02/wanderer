import sys
import random

assert (len(sys.argv) == 3)

inter_count = int(sys.argv[1])
alleys_count = int(sys.argv[2])

assert (alleys_count >= inter_count - 1)

print(inter_count, alleys_count)

alleys = []
max_alley_len = 10
repeat = False
idx = 1

for _ in range(1, alleys_count+1):
    if repeat:
        idx = random.randint(1, inter_count)
    else:
        if idx > inter_count:
            repeat = True
            idx = random.randint(1, inter_count)

    rand_inter = random.randint(1, inter_count)
    while idx == rand_inter or (idx, rand_inter) in alleys or (rand_inter, idx) in alleys:
        rand_inter = random.randint(1, inter_count)

    alleys.append((idx, rand_inter))

    rand_len = random.randint(1, max_alley_len)
    print(idx, rand_inter, rand_len)
    idx += 1

wells = []
max_wells = (inter_count // 4)
wells_count = random.randint(1, max_wells)
for i in range(1, wells_count+1):
    rand_well = random.randint(1, inter_count)
    while rand_well in wells:
        rand_well = random.randint(1, inter_count)

    wells.append(rand_well)

print()
print(wells_count, end=" ")
for well in wells:
    print(well, end=" ")
print()

exits = []
max_exits = (inter_count // 4)
exit_count = random.randint(1, max_exits)
for i in range(1, exit_count+1):
    rand_exit = random.randint(1, inter_count)
    while rand_exit in exits or rand_exit in wells:
        rand_exit = random.randint(1, inter_count)

    exits.append(rand_exit)

print(exit_count, end=" ")
for ex in exits:
    print(ex, end=" ")
print()

start = random.randint(1, inter_count)
while start in wells or start in exits:
    start = random.randint(1, inter_count)

print(1, start)

trashcans = []
max_trashcan = inter_count // 2
trashcan_count = random.randint(0, max_trashcan)
for i in range(1, trashcan_count+1):
    rand_trashcan = random.randint(1, inter_count)
    while rand_trashcan in trashcans:
        rand_trashcan = random.randint(1, inter_count)

    trashcans.append(rand_trashcan)

print(trashcan_count, end=" ")
for tr in trashcans:
    print(tr, end=" ")
print()
