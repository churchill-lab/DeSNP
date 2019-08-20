# -*- coding: utf-8 -*-

import numpy as np


class MedPolish:
    RETSTR = "overall={}\nrow={}\ncol={}\nresiduals=\n{}"

    def __init__(self, overall, row, col, residuals):
        self.overall = overall
        self.row = row
        self.col = col
        self.residuals = residuals
    
    def __str__(self):
        return self.RETSTR.format(self.overall, self.row, self.col,
                                  self.residuals)


def medpolish(x, eps=0.01, maxiter=10):
    nr, nc = x.shape
    t = 0
    r = np.repeat(0, nr)
    c = np.repeat(0, nc)
    oldsum = 0

    for iter in range(0, maxiter):
        rdelta = np.median(x, 1)
        x = x - rdelta.reshape(-1, 1)
        r = r + rdelta
        delta = np.median(c)
        c = c - delta
        t = t + delta
        cdelta = np.median(x, 0)
        x = x - cdelta
        c = c + cdelta
        delta = np.median(r)
        r = r - delta
        t = t + delta
        newsum = np.sum(np.abs(x))

        converged = newsum == 0 or abs(newsum - oldsum) < eps * newsum

        if converged:
            break
        else:
            oldsum = newsum
    
    return MedPolish(t, r, c, x)


def adjusted_medpolish(x, esp=0.01, maxiter=10):
    mp = medpolish(x, esp, maxiter)
    nu_col = mp.col

    for i in range(0, len(nu_col)):
        nu_col[i] = nu_col[i] + mp.overall

    mp.col = nu_col

    return mp


