import random
import sys
import textwrap

def main():
    if invalidArguments():
        return

    dnaDict = genDNASequence(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))

    writeSequencesToFile(dnaDict, "output.fasta")

def invalidArguments():
    if len(sys.argv) != 4:
        print("Invalid argument length: <sequenceNumber> <sequenceLengthStart> <sequenceLengthEnd>")
        return True

    try: 
        int(sys.argv[1])
        int(sys.argv[2])
        int(sys.argv[3])
    except ValueError:
        print("Invalid argument type: <int> <int> <int>")
        return True

    return False

def writeSequencesToFile(dictionary, filePath):
    with open (filePath, "w") as f:
        for k, v in dictionary.items():
            f.write(k)
            f.write('\n')
            f.write('\n'.join(textwrap.wrap(v, 80)))
            f.write('\n')

def genDNASequence(sequences, lengthStart, lengthEnd):
    dnaDict = { }
    for i in range(sequences):
        dnaSequence = ""
        for _ in range(random.randint(lengthStart, lengthEnd)):
            dnaSequence += genDNA(random.randint(0,3))
        dnaTitle = ">Sequence " + str(i + 1) + " of " + str(sequences)
        dnaDict[dnaTitle] = dnaSequence

    return dnaDict

def genDNA(x):
    return {
        0: 'A',
        1: 'C',
        2: 'G',
        3: 'T'
    }[x]

if __name__ == "__main__":
    main()
