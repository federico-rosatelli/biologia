from matplotlib import pyplot as plt
import pandas as pd


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    WARN_BOX = WARNING + '[!] '
    OK_BOX = OKBLUE + '[*] '


def printTable(table,gene,trace=[]):
    #print([i for i in range(len(table[0]))])
    indx = "      "
    for i in range(len(table[0])-1):
        indx += gene[1][i]+" "*3
    print(indx)
    print('_'*len(table[0])*4)
    for i in range(len(table)):
        print("| ",end="")
        for j in range(len(table[i])):
            #print(table[i][j])
            if (i,j) in trace:
                print(bcolors.OKGREEN+str(table[i][j])+bcolors.ENDC+" "*(2-len(str(table[i][j]))),end="| ")
            else:
                print(str(table[i][j])+" "*(2-len(str(table[i][j]))),end="| ")
        if i != 0:
            print(gene[0][i-1])
        else:
            print("")

def saveTable(table,trace=[]):
    color = [["w" for j in range(len(table[0]))]for i in range(len(table))]
    
    for point in trace:
        color[point[0]][point[1]] = "#56b5fd"
    
    fig,ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame(table)
    ax.table(cellText=df.values,cellColours=color,colLabels=df.columns,loc='center')
    #fig.tight_layout()
    plt.savefig('table.png',bbox_inches='tight')
    #plt.show()

def local_align(seq1, seq2, gap=-1,show=False)->tuple:
    m = len(seq1)
    n = len(seq2)
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    max_score = 0
    max_index = (0,0)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score = 1 if seq1[i-1] == seq2[j-1] else -1

            score_matrix[i][j] = max(
                0,
                score_matrix[i-1][j-1] + score,
                score_matrix[i-1][j] + gap,
                score_matrix[i][j-1] + gap
            )

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_index = (i, j)

    aligned_seq1 = ""
    aligned_seq2 = ""

    i, j = max_index
    traceback = []
    print(f"Table length. First sequence:{len(seq1)}, Second sequence:{len(seq2)}")
    print(f"Max term: {score_matrix[i][j]}; in index: {max_index}")
    while score_matrix[i][j] != 0:
        traceback.append((i,j))
        # if score_matrix[i-1][j-1] == 0:
        #     aligned_seq1 = seq1[i-1] + aligned_seq1
        #     aligned_seq2 = seq2[j-1] + aligned_seq2
        #     traceback.append((i-1,j-1))
        #     break
        if score_matrix[i-1][j-1] >= score_matrix[i-1][j] and score_matrix[i-1][j-1] >= score_matrix[i][j-1]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i, j = i-1, j-1
        elif score_matrix[i-1][j] >= score_matrix[i-1][j-1] and score_matrix[i-1][j] >= score_matrix[i][j-1]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    if show and len(seq1) < 80 and len(seq2) < 40:
        printTable(score_matrix,(seq1,seq2),trace=traceback)

    return aligned_seq1, aligned_seq2