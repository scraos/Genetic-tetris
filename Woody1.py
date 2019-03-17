import numpy as np
import random as rnd
import time
from copy import deepcopy
import ast
import matplotlib.pyplot as plt
from tkinter import *
from itertools import permutations


width = 10
w1 = 30
h1 = 30
grid = [[0 for i in range(0, width)] for i in range(0, width)]


def show_grid(gridval):
    for i in gridval:
        print(i)
    print()


figures = [[[1]], [[1, 1]], [[1], [1]], [[1], [1], [1]], [[1, 1], [0, 1]], [[0, 1], [1, 1]], [[1, 1], [1, 0]],
           [[1, 1, 1, 1]], [[1, 0], [1, 1], [0, 1]], [[1], [1], [1], [1]], [[1, 0], [1, 1]], [[1, 1, 1]],
           [[1,1,1],[1,0,0],[1,0,0]],[[0,0,1],[0,0,1],[1,1,1]],[[0,1],[1,1],[1,0]], [[1,1],[0,1],[0,1],[0,1]],
           [[1, 1], [1, 0], [1, 0], [1, 0]],[[1, 0], [1, 0], [1, 0], [1, 1]],[[0, 1], [0, 1], [0, 1], [1, 1]],
           [[1,1,1],[1,1,1],[1,1,1]],[[1,0,0],[1,0,0],[1,1,1]],[[1,1,1],[0,0,1],[0,0,1]],[[0,1,1],[1,1,0]],
           [[1,1,0],[0,1,1]],[[0,1,0],[1,1,1]],[[1,1,1],[0,1,0]],[[1,0],[1,1],[1,0]],[[0,1],[1,1],[0,1]],
           [[1,1],[0,1],[0,1]],[[1,1],[1,0],[1,0]],[[0,1],[0,1],[1,1]],[[1,0],[1,0],[1,1]],[[1,1],[1,1]],
           [[1,1,1,1,1]],[[1],[1],[1],[1],[1]]]
# print(len(figures))

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
            holes = countholes(curgrid)
            neighbs = count_nulneighbs(curgrid)
            if bool_place_fig(fig, [i,j], curgrid):
                place_fig(fig, [i,j], curgrid)
                totpen += a * clear_grid(curgrid) * np.sum(np.array(curgrid)) ** b
                totpen += c * (countholes(curgrid) - holes)
                totpen += d * (count_nulneighbs(curgrid)-neighbs)
                if totpen > maxscore:
                    maxscore = totpen
                    maxpos = [i,j]
    return maxpos


def genfig():
    return rnd.randint(0,len(figures)-1)


def genpop(size):
    pop = []
    for k in range(size):
        pop.append([5+rnd.random()*10,2+rnd.random()*5,-2+rnd.random()*-10,-2+rnd.random()*-8])
    return pop


def genomeplay(genome, figs):
    yes = True
    playgrid = deepcopy(grid)
    totscore = 0
    steps = 0
    while yes:
        curfig = figs[steps]
        if not check_if_pos(curfig,playgrid):
            yes = False
        toppos = genmaxscore(playgrid,curfig,genome[0],genome[1],genome[2],genome[3])
        place_fig(curfig,toppos,playgrid)
        totscore += clear_grid(playgrid)
        # show_grid(playgrid)
        # visgrid(playgrid,steps)
        steps += 1
    return steps


# def genmaxscore2(gridt, figs3, a, b, c, d):
#     maxscore = -10e6
#     maxpos = [0,0]
#     perms = list(permutations(figs3))
#     for h in perms:
#         for i in range(collen):
#             for j in range(rowlen):
#                 for k in range(collen):
#                     for l in range(rowlen):
#                         for m in range(collen):
#                             for n in range(rowlen):
#                                 totpen = 0
#                                 curgrid = deepcopy(gridt)
#                                 print(i)
#                                 if bool_place_fig(h[0], [i,j], curgrid):
#                                     if bool_place_fig(h[1], [k,l], curgrid):
#                                         if bool_place_fig(h[2], [m,n], curgrid):
#                                             totpen += a * clear_grid(curgrid) * np.sum(np.array(curgrid)) ** b
#                                             totpen += c * countholes(curgrid)
#                                             totpen += d * count_nulneighbs(curgrid)
#                                     if totpen > maxscore:
#                                         maxscore = totpen
#                                         maxpos = [[i,j],[k,l],[m,n]]
#     return maxpos
#
#
# def genomeplay2(genome, figslist):
#     yes = True
#     playgrid = deepcopy(grid)
#     # totscore = 0
#     steps = 0
#     while yes:
#         curfig = figslist[steps]
#         if not check_if_pos(curfig[0],playgrid):
#             yes = False
#         elif not check_if_pos(curfig[1],playgrid):
#             yes = False
#         elif not check_if_pos(curfig[2],playgrid):
#             yes = False
#         toppos = genmaxscore2(playgrid,curfig,genome[0],genome[1],genome[2],genome[3])
#         place_fig(curfig[0], toppos[0], playgrid)
#         place_fig(curfig[1], toppos[1], playgrid)
#         place_fig(curfig[2], toppos[2], playgrid)
#         # totscore += clear_grid(playgrid)
#         # show_grid(playgrid)
#         visgrid(playgrid,steps)
#         steps += 1
#     return steps

def showpams(pams,avgs,t,stdd):
    plt.subplot(3, 2, 1)
    plt.hist(pams[0],bins=30)
    plt.subplot(3, 2, 2)
    plt.hist(pams[1],bins=30)
    plt.subplot(3, 2, 3)
    plt.hist(pams[2],bins=30)
    plt.subplot(3, 2, 4)
    plt.hist(pams[3],bins=30)
    plt.subplot(3,2,5)
    plt.plot(np.linspace(1,t,t),avgs)
    plt.plot(np.linspace(1, t, t),stdd)
    plt.show()

def appenpams(pam,pamnum,pams):
    for i in pams:
        pam.append(i[pamnum])


def mutpop(pop,top):
    time = 1
    newpop = pop
    a = []
    b = []
    c = []
    d = []
    avg = []
    stddevs = []
    mutrate = 0.2
    while time < 50:
        oldpop = deepcopy(newpop)
        results = []
        figs = [rnd.randint(0,len(figures)-1) for i in range(1000)]
        for k in oldpop:
            results.append(genomeplay(k,figs))
            # print(k)
        topn = np.argpartition(-np.array(results), top)[:top]
        print(time)
        print(results)
        print(topn)
        print('ensemble avg {}'.format(np.mean(results)))
        print('ensemble stddev {}'.format(np.std(results)))
        print('lowest performance')
        print(min(results))
        print('highest performance')
        print(max(results))
        print()
        avg.append(np.mean(np.array(results)))
        stddevs.append(np.std(results))
        appenpams(a, 0, oldpop)
        appenpams(b, 1, oldpop)
        appenpams(c, 2, oldpop)
        appenpams(d, 3, oldpop)
        if time % 5 == 0:
            showpams([a,b,c,d],avg,time,stddevs)
            mutrate -= (time/50)*mutrate
        newpop = []
        # better genomes should have relatively more offspring
        # past performance should be taken into account, if avg and stddev are both low
        for j in topn:
            newpop.append([oldpop[j][0], oldpop[j][1], oldpop[j][2], oldpop[j][3]])
            for m in range(len(pop)//top-1):
                newpop.append([oldpop[j][0]*rnd.gauss(1,mutrate),oldpop[j][1]*rnd.gauss(1,mutrate),
                               oldpop[j][2]*rnd.gauss(1,mutrate),oldpop[j][3]*rnd.gauss(1,mutrate)])
        time += 1
    return pop

mutpop(genpop(50),25)
# root = Tk()
# w = Canvas(root, width=500, height=400)
# w.pack()
#
# steps = Label(root, text='0')
# steps.pack(side=RIGHT)
#
# notation = [None for i in range(20)]
# for i in range(10):
#     notation[i] = Label(root,text='{}'.format(i))
#     notation[i].place(x=3+i*w1,y = width*h1+2)
#
# for i in range(10):
#     notation[i+10] = Label(root, text='{}'.format(i))
#     notation[i+10].place(x=width*w1+2,y=2+i*h1)
#
# def visgrid(grids,step):
#     w.delete("all")
#     time.sleep(0.05)
#     steps.config(text='{}'.format(step))
#     for i in range(10):
#         for j in range(10):
#             if grids[i][j] == 0:
#                 w.create_rectangle(j*w1, i*h1, (j+1)*w1, (i+1)*h1, fill='white')
#             if grids[i][j] == 1:
#                 w.create_rectangle(j * w1, i * h1, (j+1) * w1, (i+1) * h1, fill='#404040')
#     root.update()
#
# pams = [11,2,-3,-2]
# enterfig = Entry(root)
# enterfig.pack(side = RIGHT)
# figtext = Label(root, text='Calc for fig:')
# figtext.pack(side = RIGHT)
# playgrid1 = deepcopy(grid)
# bestmovetext = Label(root, text='')
# var1 = IntVar()
# def prbmove():
#     bestmovetext.config(text='{}'.format(genmaxscore(playgrid1,int(enterfig.get()),pams[0],pams[1],
#                                                                         pams[2],pams[3])))
#     var1.set(1)
#
#
# printbestmove = Button(root, text='Print bestmove', command=prbmove)
# printbestmove.pack(side=LEFT)
# bestmovetext.pack(side = LEFT)
# root.update()
# stp = 0
# while True:
#     visgrid(playgrid1, stp)
#     printbestmove.wait_variable(var1)
#     b = ast.literal_eval(bestmovetext.cget("text"))
#     place_fig(int(enterfig.get()),b,playgrid1)
#     clear_grid(playgrid1)
#     visgrid(playgrid1,stp)
#     stp += 1

# def drawfig(fig,pos):
#     for i in range(give_row(fig)):
#         for j in range(give_col(fig)):
#             if figures[fig][i][j] == 0:
#                 w.create_rectangle(pos[1]*w1+j*w1, pos[0]*h1+i*h1, pos[1]*w1+(j+1)*w1,
#                                    pos[0]*h1+(i+1)*h1, fill='white')
#             if figures[fig][i][j] == 1:
#                 w.create_rectangle(pos[1]*w1+j*w1, pos[0]*h1+i*h1, pos[1]*w1+(j+1)*w1,
#                                    pos[0]*h1+(i+1)*h1, fill='#404040')
#
# fignums = [None for i in range(len(figures))]
# for i in  range(7):
#     drawfig(i,[0,2+6*i])
#     fignums[i] = Label(root,text=i)
#     fignums[i].place(y=0,x=(6*i)*h1)
# for i in  range(7):
#     drawfig(i+7,[6,2+6*i])
#     fignums[i+7] = Label(root, text=i+7)
#     fignums[i+7].place(y=6*w1, x=(6*i)*h1)
# for i in  range(7):
#     drawfig(i+14,[12,2+6*i])
#     fignums[i+14] = Label(root, text=i+14)
#     fignums[i+14].place(y=12*w1,x=(6*i)*h1)
# for i in  range(7):
#     drawfig(i+21,[18,2+6*i])
#     fignums[i+21] = Label(root, text=i+21)
#     fignums[i+21].place(y=18*w1, x=(6*i)*h1)
# for i in  range(7):
#     drawfig(i+28,[24,2+6*i])
#     fignums[i+28] = Label(root, text=i+28)
#     fignums[i+28].place(y=24*w1, x=(6*i)*h1)
# root.mainloop()

# figs = [rnd.randint(0,len(figures)-1) for i in range(1000)]
# genomeplay(pams,figs)