#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, July 2014.

from Tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import sys
import random

class Percolation:

    def __init__(self, L=61):
        if L % 2 == 0:
            raise ValueError("lattice size L must be odd number.")
        self.sub = None
        self.lattice = None
        self.L = L # lattice size

    def perc_cluster(self):
        self.lattice = np.zeros([self.L+2, self.L+2], dtype=int)
        self.lattice[:1,:] = self.lattice[:, :1] = -1
        self.lattice[self.L+1:,:] = self.lattice[:, self.L+1:] = -1
        center = (self.L//2) + 1
        self.lattice[center, center] = 1
        nextseed = [(center, center)]
        if self.sub is None or not self.sub.winfo_exists():
            lattice = self.lattice
            L = self.L
            rn = np.random.random
            choice = random.choice
            ne = [(0, -1), (0, 1), (-1, 0), (1, 0)]
            nnsite = set([(center+nx, center+ny) for nx, ny in ne])
            t = [0] # time
            S = [4] # a number of growing sites
            N = [1] # a number of occupied sites
            percolate = False
            l = set([])
            while percolate == False:
                nn = choice(list(nnsite))
                nnsite.remove(nn)
                lattice[nn] = 1
                i, j = nn
                nnsite = nnsite | set([(i+nx, j+ny) for nx, ny in ne
                                        if lattice[i+nx, j+ny] == 0])
                if i == 1:
                    l.add('top')
                if i == L:
                    l.add('bottom')
                if j == 1:
                    l.add('left')
                if j == L:
                    l.add('right')
                
                glsite = (np.array([a[0] for a in list(nnsite)]), 
                          np.array([a[1] for a in list(nnsite)]))
                lattice[glsite] = -1
                if ('top' in l and 'bottom' in l) or \
                   ('right' in l and 'left' in l):
                    percolate = True
                
                t.append(t[-1]+1)
                S.append(len(nnsite))
                N.append(np.sum(lattice == 1))
            self.lattice = lattice[1:-1, 1:-1]
        
        return t, S, N

    def draw_canvas(self, rect, L, show=1):
        default_size = 640 # default size of canvas
        r = int(default_size/(2*L))
        if r == 0:
            r = 1
        fig_size = 2*r*L
        margin = 10
        sub = Toplevel()

        self.canvas = Canvas(sub, width=fig_size+2*margin,
                             height=fig_size+2*margin)
        self.canvas.create_rectangle(margin, margin,
                                    fig_size+margin, fig_size+margin,
                                    outline='black', fill='white')
        self.canvas.pack()

        c = self.canvas.create_rectangle

        site = np.where(rect == show) # 1: occupied site, -1: growing site
        for m, n in zip(site[0], site[1]):
            c(2*m*r+margin, 2*n*r+margin,
              2*(m+1)*r+margin, 2*(n+1)*r+margin,
              outline='', fill='black')

class TopWindow:

    def quit(self):
        self.root.destroy()
        sys.exit()

    def show_window(self, pr, pushed, b4_pushed, b5_pushed):
        self.root = Tk()
        self.root.title('Percolation')
        f1 = Frame(self.root, padx=5, pady=5)
        b1 = Button(f1, text='run (and show figure)', command=pushed)
        b1.pack(side='top', expand=YES, fill='x')
        f1.pack(fill='x')
        
        f2 = Frame(self.root, padx=5, pady=5)
        b5 = Button(f2, text='show growing site', command=b5_pushed)
        b5.pack(side='top', expand=YES, fill='x')

        b4 = Button(f2, text='plot graph', command=b4_pushed)
        b4.pack(expand=YES, fill='x')

        b2 = Button(f2, text='save canvas to sample.eps', command=pr)
        b2.pack(expand=YES, fill='x')
        f2.pack(fill='x')

        f3 = Frame(self.root, padx=5, pady=5)
        b3 = Button(f3, text='quit', command=self.quit)
        b3.pack(expand=YES, fill='x')
        f3.pack(fill='x')

        self.root.mainloop()


if __name__ == '__main__':
    L = 201
    top = TopWindow()
    per = Percolation(L=L)
    count = 1

    def pr():
        global count
        d = per.canvas.postscript(file="figure_%d.eps" % count)
        print "saved the figure to a eps file"
        count += 1

    def pushed():
        global t, S, N
        t, S, N = per.perc_cluster()
        per.draw_canvas(per.lattice, L)

    def b5_pushed():
        if per.lattice == None:
            per.perc_cluster()
        per.draw_canvas(per.lattice, L, show=-1)

    def b4_pushed():
        if per.lattice == None:
            per.perc_cluster()
        sqrt = np.sqrt
        c = (L//2)
        rlattice = np.array([sqrt(x**2 + y**2) for x in range(c) 
                            for y in range(c)]).reshape(c, c)
        i = 1
        _r = 2
        r, M = [], []
        while _r <= c:
            s = np.where(rlattice <= _r)
            set_in_r = set(zip(np.append(s[0], [s[0], -s[0], -s[0]])+c,
                               np.append(s[1], [-s[1], s[1], -s[1]])+c))
            l = list(set_in_r)
            site = (np.array([p[0] for p in l]), np.array([p[1] for p in l]))
            M_r = np.sum(per.lattice[site] == 1)
            r.append(_r)
            M.append(M_r)
            i += 1
            _r = 2**i
            
        log = np.log
        def fit_func(parameter0, r, M_r):
            c1 = parameter0[0]
            c2 = parameter0[1]
            residual = log(M_r) - c1 - c2*log(r)
            return residual
        
        r = np.array(r)
        M = np.array(M)
        x = np.logspace(np.log10(np.min(r))-0.1, np.log10(np.max(r))+0.1, 10)
        parameter0 = [0.1, 2.0]
        result = optimize.leastsq(fit_func, parameter0, args=(r, M))
        c1 = result[0][0]
        D = result[0][1]

        def fitted(r, c1, D):
            return np.exp(c1)*(r**D)
        
        fx = fitted(x, c1, D)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.14)
        plt.subplots_adjust(right=0.68)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$\ln r$', fontsize=16)
        ax.set_ylabel(r'$\ln M$', fontsize=16)
        ax.set_ymargin(0.05)
        ax.plot(r, M, 'o')
        ax.plot(x, fx, '-', color='black',
                label='\n$\mathrm{fitted\ curve}$\n$(D= %1.2f)$' % D)
        ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', borderaxespad=0)
        plt.show()
        
    top.show_window(pr, pushed, b4_pushed, b5_pushed)

