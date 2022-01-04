from Bio import SeqIO
from Bio.Seq import Seq
from numpy import around
start_codons = ['ATG']
stop_codons = ['TAG', 'TAA', 'TGA']


def read_fasta(handle: str):
    for seq_record in SeqIO.parse(handle, "fasta"):
        return seq_record


def get_frame(seq: Seq, i: int):
    return seq[i::3]


def find_start_to_stop(frame, stride, x, length):
    found = False
    for j in range(len(frame[stride:])):
        x += frame[j + stride]
        length += 1
        if frame[j + stride:j + stride + 3] in stop_codons:
            x += frame[j + stride + 1] + frame[j + stride + 2]
            found = True
            length += 2
            break
    return found, x


def find_stop_to_start(frame, stride, i, length):
    found = False
    x = ""
    for j in range(len(frame[stride:])):
        if frame[j + stride: j + stride + 3] in start_codons:
            x = frame[i: j + stride + 3]
            length += 3
            found = True
            j += 3
        elif frame[j + stride: j + stride + 3] in stop_codons:
            break
        length += 1
    return found, x


def find_codons(arr: list, frame: Seq, inArr: list, startOrStop: bool):
    for i in range(len(frame)):
        stride = i + 3
        length = 0
        if frame[i:stride] in inArr:
            length += 3
            if startOrStop:
                found, x = find_start_to_stop(
                    frame, stride, frame[i:stride], length)
            else:
                found, x = find_stop_to_start(frame, stride, i, length)
            if found:
                arr.append(x)
        i += length


def find_pairs(seq: Seq, startOrStop):
    arr = []
    for i in range(3):
        frame = get_frame(seq, i)
        if startOrStop:
            find_codons(arr, frame, start_codons, startOrStop)
        else:
            find_codons(arr, frame, stop_codons, startOrStop)
    rcSeq = seq.reverse_complement()
    for i in range(3):
        frame = get_frame(rcSeq, i)
        if startOrStop:
            find_codons(arr, frame, start_codons, startOrStop)
        else:
            find_codons(arr, frame, stop_codons, startOrStop)
    return arr


def filter_length(arr: list):
    return [b for b in arr if len(b) >= 100]


def get_possible(seq: Seq, n1: int, n2: int):
    arr = []
    i = 0
    while i < len(seq):
        if seq[i:i+n2] not in arr:
            if len(seq[i:i+n2]) == n2:
                arr.append(seq[i:i+n2])
        i += n1
    return arr


def concat_sequences(seqs):
    sum = Seq("")
    for item in seqs:
        sum += item
    return sum


def get_freqs(seq: Seq, n1: int, n2: int):
    freqs = []
    allPossible = get_possible(seq, n1, n2)
    for possible in allPossible:
        freqs.append(Frequency(possible, seq.count(possible) / len(seq)))
    return freqs


def get_genomes(genomes: int):
    arr = []
    type = "bacterial"
    for i in range(genomes):
        num = i + 1
        if num > genomes / 2:
            type = "mamalian"
            num -= int(genomes / 2)
        arr.append(read_fasta("data/" + type + str(num) + ".fasta"))
    return arr


def get_pairs(genomes, startOrStop):
    arr = []
    for i in range(len(genomes)):
        arr.append(find_pairs(genomes[i].seq, startOrStop))
    return arr


class Frequency:
    def __init__(self, code, freq) -> None:
        self.code = code
        self.freq = freq


def get_freqs_arr(genomes, codes1, codes2, n1, n2):
    arr = []
    for i in range(len(genomes)):
        concated = concat_sequences(codes1[i]) + concat_sequences(codes2[i])
        if concated:
            arr.append(get_freqs(concated, n1, n2))
        else:
            arr.append([])
    return arr


def find_same_freq(code, arr: list):
    for x in arr:
        if x.code == code:
            return x


def compare_freqs(arr1: list, arr2: list):
    resultArr = []
    for i in range(len(arr1)):
        sameCodeItem = find_same_freq(arr1[i].code, arr2)
        if sameCodeItem:
            delta = arr1[i].freq - sameCodeItem.freq
            if delta < 0:
                delta = delta * -1
            resultArr.append(delta)
    if len(resultArr) > 0:
        return (sum(resultArr) / len(resultArr)) * 1000
    return -1


def create_philyp_matrix(genomes, freqs):
    philyp_matrix = []
    for i in range(len(genomes)):
        philyp_matrix.append([])
    for i in range(len(genomes)):
        for j in range(len(genomes)):
            philyp_matrix[i].append(compare_freqs(
                freqs[i], freqs[j]))
    return philyp_matrix


def print_philyp_matrix(genomes, matrix):
    print(len(genomes))
    for i in range(len(genomes)):
        newList = around(matrix[i], 4)
        print(genomes[i].id, *newList)


# setup

genomes = get_genomes(8)

# 1

startToStop = get_pairs(genomes, True)

# 2

stopToStart = get_pairs(genomes, False)

# 3

for i in range(len(genomes)):
    startToStop[i] = filter_length(startToStop[i])
    stopToStart[i] = filter_length(stopToStart[i])

# 4

codonFreqsArr = get_freqs_arr(genomes, startToStop, stopToStart, 1, 3)
dicodonFreqsArr = get_freqs_arr(genomes, startToStop, stopToStart, 3, 6)

# 5

print_philyp_matrix(genomes, create_philyp_matrix(genomes, codonFreqsArr))
print_philyp_matrix(genomes, create_philyp_matrix(genomes, dicodonFreqsArr))
