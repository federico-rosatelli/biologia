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
    indx = "  "
    for i in range(trace[len(trace)-1][1],trace[len(trace)-1][1]+10):
        indx += gene[1][i]+" "*3
    print(indx)
    print('_'*10*4)
    for i in range(trace[len(trace)-1][0],trace[len(trace)-1][0]+10):
        print("| ",end="")
        for j in range(trace[len(trace)-1][1],trace[len(trace)-1][1]+10):
            #print(table[i][j])
            if (i,j) in trace:
                print(bcolors.OKGREEN+str(table[i][j])+bcolors.ENDC+" "*(2-len(str(table[i][j]))),end="| ")
            else:
                print(str(table[i][j])+" "*(2-len(str(table[i][j]))),end="| ")
        if i != 0:
            print(gene[0][i-1])
        else:
            print("")
    print('\n')

def saveTable(table,trace=[]):
    tableCopy = []
    for i in range(10,40):
        tbCp = []
        for j in range(10,40):
            tbCp.append(table[i][j])
        tableCopy.append(tbCp)
    color = [["w" for j in range(len(tableCopy[0]))]for i in range(len(tableCopy))]
    for i in range(len(tableCopy)):
        for j in range(len(tableCopy[i])):
            if (i+10,j+10) in trace:
                color[i][j] = "#56b5fd"
    # for i in range(10,20):
    #     color[trace[len(trace)-i-1][0]][trace[len(trace)-i-1][1]] = "#56b5fd"
    
    fig,ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame(tableCopy)
    ax.table(cellText=df.values,cellColours=color,colLabels=df.columns,loc='center')
    #fig.tight_layout()
    plt.savefig('table.png',bbox_inches='tight')
    #plt.show()


class Alignment:
    def __init__(self,*seqs:str,gap:int=1,show_table:bool=False) -> None:
        if len(seqs)<2:
            print(bcolors.FAIL+f"MinMaxError: input len for sequences must be >= 2, got {len(seqs)}"+bcolors.ENDC)
            return
        self.seqs = seqs
        self.seq1 = seqs[0]
        self.seq2 = seqs[1]
        self.gap = gap
        self.show_table = show_table
        self.aligned_seq1 = str
        self.aligned_seq2 = str
    
    def __str__(self) -> str:
        return self.seq1 + "\n" + self.seq2

    def __len__(self) -> int:
        return len(self.seq1)*len(self.seq2)

    def createScoreMatrix(self,lnSeq1:int,lnSeq2:int) -> tuple:
        score_matrix = [[0 for _ in range(lnSeq2 + 1)] for _ in range(lnSeq1 + 1)]

        max_score = 0
        max_index = (0,0)
        for i in range(1, lnSeq1 + 1):
            for j in range(1, lnSeq2 + 1):
                score = 1 if self.seq1[i-1] == self.seq2[j-1] else -1

                score_matrix[i][j] = max(
                    0,
                    score_matrix[i-1][j-1] + score,
                    score_matrix[i-1][j] - self.gap,
                    score_matrix[i][j-1] - self.gap
                )

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_index = (i, j)
        return score_matrix,max_index

    def localAlignment(self,save_table=False) -> tuple:
        if len(self.seqs)>2:
            print(bcolors.WARN_BOX+f"Warning! Only 2 arguments were expected, but got {len(self.seqs)}.\n\t-The algorithm will use only the first 2 sequences..."+bcolors.ENDC)
        aligned_seq1 = ""
        aligned_seq2 = ""

        score_matrix, max_index_score_matrix = self.createScoreMatrix(len(self.seq1),len(self.seq2))

        i, j = max_index_score_matrix
        traceback = []
        print(f"Table length. First sequence:{len(self.seq1)}, Second sequence:{len(self.seq2)}")
        print(f"Max term: {score_matrix[i][j]}; in index: {max_index_score_matrix}")
        f1 = 0
        while score_matrix[i][j] != 0:
            traceback.append((i,j))
            if score_matrix[i-1][j-1] >= score_matrix[i-1][j] and score_matrix[i-1][j-1] >= score_matrix[i][j-1]:
                aligned_seq1 = self.seq1[i-1] + aligned_seq1
                aligned_seq2 = self.seq2[j-1] + aligned_seq2
                i, j = i-1, j-1
                f1 += 1
            elif score_matrix[i-1][j] >= score_matrix[i-1][j-1] and score_matrix[i-1][j] >= score_matrix[i][j-1]:
                aligned_seq1 = self.seq1[i-1] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2
                i -= 1
            else:
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = self.seq2[j-1] + aligned_seq2
                j -= 1
        
        if self.show_table:
            printTable(score_matrix,(self.seq1,self.seq2),trace=traceback)
            if save_table:
                saveTable(score_matrix,trace=traceback)
        print(f"{(f1/len(aligned_seq1))*100}% of alignment")
        return aligned_seq1, aligned_seq2
    
    def globalAlignment(self):
        return 0