import numpy as np
import random as rnd
from copy import deepcopy


width = 10
grid = [[0 for i in range(0, width)] for i in range(0, width)]


def show_grid(gridval):
    for i in gridval:
        print(i)
    print()


figures = [[[1]], [[1, 1]], [[1], [1]], [[1], [1], [1]], [[1, 1], [0, 1]], [[0, 1], [1, 1]], [[1, 1], [1, 0]],
           [[1, 1, 1, 1]], [[1, 0], [1, 1], [0, 1]], [[1], [1], [1], [1]], [[1, 0], [1, 1]], [[1, 1, 1]],
           [[1,1,1],[1,0,0],[1,0,0]],[[0,0,1],[0,0,1],[1,1,1]],[[1,0],[1,1],[0,1]], [[1,1],[0,1],[0,1],[0,1]],
           [[1, 1], [1, 0], [1, 0], [1, 0]],[[1, 0], [1, 0], [1, 0], [1, 1]],[[0, 1], [0, 1], [0, 1], [1, 1]],
           [[1,1,1],[1,1,1],[1,1,1]]]


def give_row(fig):
    return len(figures[fig])


def give_col(fig):
    return len(figures[fig][0])


def bool_place_fig(fig, pos, gridt):
    if (give_col(fig)+pos[1] > width) or (pos[1] < 0):
        return False
    if (give_row(fig)+pos[0] > width) or (pos[0] < 0):
        return False
    curposs = 0
    for i in range(give_row(fig)):
        for j in range(give_col(fig)):
            if figures[fig][i][j] == 1:
                curposs += gridt[pos[0] + i][pos[1] + j]
            if curposs > 0:
                return False
    return True

def place_fig(fig,pos,gridt):
    for i in range(give_row(fig)):
        for j in range(give_col(fig)):
            if figures[fig][i][j] == 1:
                gridt[pos[0] + i][pos[1] + j] = 1
    return True

collen = len(grid)
rowlen = len(grid[0])

def check_if_pos(fig, gridt):
    for i in range(collen):
        for j in range(rowlen):
            curgrid = deepcopy(gridt)
            if bool_place_fig(fig,[i,j],curgrid):
                return True
    return False

def transpose_grid(gridt):
    return list(map(list, zip(*gridt)))


def clear_grid(curgrid):
    addscore = 0
    for i in range(collen):
        if 0 not in curgrid[i]:
            curgrid[i] = [0 for j in range(0, width)]
            addscore += 16
    gridt = transpose_grid(curgrid)
    for i in range(rowlen):
        if 0 not in gridt[i]:
            for j in range(collen):
                curgrid[j][i] = 0
            addscore += 16
    return addscore


def giveneighbs(pos):
    neighbs = []
    if pos[0] > 0:
        neighbs.append([pos[0]-1,pos[1]])
    if pos[0] < collen - 1:
        neighbs.append([pos[0]+1,pos[1]])
    if pos[1] > 0:
        neighbs.append([pos[0],pos[1]-1])
    if pos[1] < rowlen - 1:
        neighbs.append([pos[0],pos[1]+1])
    return neighbs


# rowlen is numcols
# collen is numrows

def countholes(gridt):
    tocheck = []
    comps = 0
    for i in range(width):
        for j in range(width):
            if gridt[i][j] == 0:
                tocheck.append([i,j])
    while tocheck:
        cur = tocheck.pop()
        comps += 1
        queue = []
        for k in giveneighbs(cur):
            if k in tocheck:
                queue.append(k)
        while queue:
            cur2 = queue.pop()
            tocheck.remove(cur2)
            for m in giveneighbs(cur2):
                if (m in tocheck) and (m not in queue):
                    queue.append(m)
    return comps


def count_nulneighbs(gridt):
    numnuls = 0
    for i in range(collen):
        for j in range(rowlen):
            if gridt[i][j] == 1:
                for k in giveneighbs([i,j]):
                    if gridt[k[0]][k[1]] == 0:
                        numnuls += 1
    return numnuls

def genmaxscore(gridt, fig, a, b, c, d):
    maxscore = -10e6
    maxpos = [0,0]
    for i in range(collen):
        for j in range(rowlen):
            totpen = 0
            curgrid = deepcopy(gridt)
            if bool_place_fig(fig, [i,j], curgrid):
                place_fig(fig, [i,j], curgrid)
                totpen += a * clear_grid(curgrid) * np.sum(np.array(curgrid)) ** b
                totpen += c * countholes(curgrid)
                totpen += d * count_nulneighbs(curgrid)
                if totpen > maxscore:
                    maxscore = totpen
                    maxpos = [i,j]
    return maxpos


def genfig():
    return rnd.randint(0,len(figures)-1)


def genpop(size):
    pop = []
    for k in range(size):
        pop.append([rnd.random()*10,rnd.random()*10,rnd.random()*-10,rnd.random()*-10])
    return pop


def genomeplay(genome):
    yes = True
    playgrid = deepcopy(grid)
    totscore = 0
    steps = 0
    while yes:
        curfig = genfig()
        if not check_if_pos(curfig,playgrid):
            yes = False
        toppos = genmaxscore(playgrid,curfig,genome[0],genome[1],genome[2],genome[3])
        place_fig(curfig,toppos,playgrid)
        totscore += clear_grid(playgrid)
        # show_grid(playgrid)
        steps += 1
    return steps

def mutpop(pop,top):
    time = 1
    newpop = pop
    while time < 100:
        oldpop = deepcopy(newpop)
        results = []
        for k in pop:
            results.append(genomeplay(k))
            # print(k)
        topn = np.argpartition(np.array(results), -top)[-top:]
        print(results)
        print([oldpop[i] for i in topn])
        mutrate = 0.2*time
        newpop = []
        for j in topn:
            for m in range(len(pop)//top):
                newpop.append([oldpop[j][0]*rnd.gauss(1,mutrate),oldpop[j][1]*rnd.gauss(1,mutrate),oldpop[j][2]*rnd.gauss(1,mutrate),oldpop[j][3]*rnd.gauss(1,mutrate)])
        time += 1
    return pop
