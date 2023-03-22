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


def createTables(gene1,gene2):
    l1 = len(gene1)
    l2 = len(gene2)
    table = []
    for i in range(l1+1):
        
        tab = []
        for k in range(l2+1):
            if i == 0:
                tab.append(0)
            elif k == 0:
                tab.append(0)
            else:
                tab.append(0)
        table.append(tab)
    gap = []
    for i in range(l1):
        g = []
        for k in range(l2):
            
            if gene1[i] == gene2[k]:
                g.append(1)
            else:
                g.append(-1)
        gap.append(g)
    return table,gap

    

def algorithm(gene1,gene2,show=False):
    """AUTHOR: Federico Rosatelli (federico-rosatelli)\n
    Smith-Waterman algorithm for local sequence alignment
    """
    table,gap = createTables(gene1,gene2)
    
    print(f"Table length. First sequence:{len(table)}, Second sequence:{len(table[0])}")
    maxim = 0
    ind = (0,0)
    gap1 = 1
    for i in range(1,len(table)):
        for j in range(1,len(table[i])):  
            a = table[i-1][j-1] + gap[i-1][j-1]
            b = table[i-1][j] - gap1
            c = table[i][j-1] - gap1
            table[i][j] = max(a,b,c,0)
            if table[i][j] >= maxim:
                maxim = table[i][j]
                ind = (i,j)
        #gap1 += 2
    print(f"Max term: {maxim}; in index: {ind}")
    trace = traceback(table,ind)
    if len(gene1) < 40 and len(gene2) < 40 and show:
    #printTable(gap,(gene1,gene2))
    #print(gene1,gene2)
        printTable(table,(gene1,gene2),trace=trace[1])
    #saveTable(table,trace=trace[1])
    g1 = ""
    g2 = ""

    route = trace[1]
    #print(route)
    for i in range(1,len(route)):
        point = route[i]
        if gene1[point[0]] == gene2[point[1]]:
            g1 += gene1[point[0]]
            g2 += gene2[point[1]]
        else:
            next_point = route[i-1]
            if next_point[0] == point[0] and next_point[1] != point[1]:
                g1 += "-"
                g2 += gene2[point[1]]
            else:
                g1 += gene1[point[0]]
                g2 += "-"
            


    return g1[::-1],g2[::-1]

def traceback(table,start):
    """Traceback algorithm based on Dijkstra's follow-the-minimum algorithm"""
    open_dict = {}
    close_list = {}
    point = (start,[start],table[start[0]][start[1]])
    open_dict[start] = point
    while open_dict:
        diag,num_diag = (point[0][0]-1,point[0][1]-1),table[point[0][0]-1][point[0][1]-1]
        # if num_diag == 0:
        #     return (diag,open_dict[point[0]][1]+[diag],open_dict[point[0]][2]+num_diag)
        up,num_up = (point[0][0]-1,point[0][1]),table[point[0][0]-1][point[0][1]]
        # if num_up == 0:
        #     return (up,open_dict[point[0]][1]+[up],open_dict[point[0]][2]+num_up)
        left,num_left = (point[0][0],point[0][1]-1),table[point[0][0]][point[0][1]-1]
        # if num_left == 0:
        #     (left,open_dict[point[0]][1]+[left],open_dict[point[0]][2]+num_left)
        
        if (diag not in open_dict or open_dict[diag][2] <  open_dict[point[0]][2]+num_diag) and diag not in close_list:
            open_dict[diag] = (diag,open_dict[point[0]][1]+[diag],open_dict[point[0]][2]+num_diag)
        if (up not in open_dict or open_dict[up][2] <  open_dict[point[0]][2]+num_up) and up not in close_list:
            open_dict[up] = (up,open_dict[point[0]][1]+[up],open_dict[point[0]][2]+num_up)
        if (left not in open_dict or open_dict[left][2] <  open_dict[point[0]][2]+num_left) and left not in close_list:
            open_dict[left] = (left,open_dict[point[0]][1]+[left],open_dict[point[0]][2]+num_left)
        
        
        
        
        close_list[point[0]] = point
        open_dict.pop(point[0],None)
        max_next_point = max(open_dict.keys(),key=lambda k: open_dict[k][2])
        #print(max_next_point)
        if max_next_point[0] == 0 or max_next_point[1] == 0:
            return open_dict[max_next_point]
        point = open_dict[max_next_point]
        
