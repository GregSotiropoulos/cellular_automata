# -*- coding: utf-8 -*-

"""
.. module:: conway_patterns
   :platform: Unix, Windows
   :synopsis: Various patterns for Conway's Game of Life

.. moduleauthor:: Greg Sotiropoulos <greg.sotiropoulos@gmail.com>

Some patterns for the Game of Life. Some well-known patterns are listed here,
such as the hallmark tiny glider and "Gosper's glider gun". This is a sort of
experimental area, and many patterns have come and gone. The result is that I
no longer remember what each and every pattern is/does. They are left here for
the curious; several of them will be of interest but others will not.
"""

import numpy as np

ui8, npr = np.uint8, np.random


# some helper functions

def npa(x):
    return np.array(x, dtype=ui8)


def npz(x):
    return np.zeros(x, dtype=ui8)


def npo(x):
    return np.ones(x, dtype=ui8)


# well-known patterns

def glider():
    """The classic glider"""
    return npa([[0, 0, 1],
                [1, 0, 1],
                [0, 1, 1]])


def gosper_glider_gun():
    """creates a Gosper Glider Gun"""
    gun = npz((11, 38))
    gun[5:7, 1:3] = gun[3:5, 35:37] = 1

    gun[3][13] = gun[3][14] = 1
    gun[4][12] = gun[4][16] = 1
    gun[5][11] = gun[5][17] = 1
    gun[6][11] = gun[6][15] = gun[6][17] = gun[6][18] = 1
    gun[7][11] = gun[7][17] = 1
    gun[8][12] = gun[8][16] = 1
    gun[9][13] = gun[9][14] = 1

    gun[1][25] = 1
    gun[2][23] = gun[2][25] = 1
    gun[3][21] = gun[3][22] = 1
    gun[4][21] = gun[4][22] = 1
    gun[5][21] = gun[5][22] = 1
    gun[6][23] = gun[6][25] = 1
    gun[7][25] = 1

    return gun


def still_elipse():
    """
    Smallest possible elliptical pattern. Turns up all the time towards
    the end of a simulation.
    """
    return npa([
        [0, 1, 1, 0],
        [1, 0, 0, 1],
        [0, 1, 1, 0]
    ])


def three_bar():
    return npa([[1, 1, 1]])


# ------- custom, tests etc -------------

def still_elipses_2_touch():
    el = still_elipse()
    el = np.pad(el, 1)
    zs = npz((5, 15))
    zs_ver = npz((5, 3))
    pat = np.concatenate((el, zs_ver, el), axis=1)
    pat = np.concatenate((zs, pat, zs), axis=0)
    pat += np.rot90(pat)
    pat = np.concatenate((pat[:, :-1], pat[:, 1:]), axis=1)
    return pat


def pi_shuttle():
    return npa([
        [1, 1, 1],
        [1, 0, 1],
        [1, 0, 1,],

    ])


def hollow_sq3_cut_corners():
    return npa([
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1]
    ])


def new_still_2():
    hsq3cc = hollow_sq3_cut_corners()
    return np.concatenate((hsq3cc,
                           npz((1, 3)),
                           np.flipud(hsq3cc)), axis=0)


def new_still_1():
    return npa([
        [0, 0, 1, 1],
        [0, 0, 1, 0],
        [1, 0, 1, 0],
        [1, 1, 0, 0]
    ])


def new_still_3():
    return npa([
        [1, 1, 0, 0],
        [1, 0, 0, 1],
        [0, 1, 1, 1],
        [0, 0, 0, 0],
        [0, 1, 1, 1],
        [1, 0, 0, 1],
        [1, 1, 0, 0]
    ])


def new_still_4():
    return npa([
        [1, 1, 0, 1],
        [1, 0, 1, 1],
    ])


def new_oscillator2_1():
    return npa([
        [0, 0, 1, 1],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [1, 1, 0, 0]
    ])


def new_oscillator15_1():
    return npa([
        [0, 1, 0, 0, 1, 0],
        [1, 0, 1, 1, 0, 1],
        [0, 1, 0, 0, 1, 0],
        [0, 1, 0, 0, 1, 0],
        [1, 0, 1, 1, 0, 1],
        [0, 1, 0, 0, 1, 0]
    ])


def osc_15_group():
    pat = new_oscillator15_1()
    return np.concatenate((pat, npz((6, 2)), pat), axis=1)


def new_oscillator6_1():
    pat = npa([
        [0, 0, 1, 0, 0],
        [1, 1, 0, 1, 1],
        [0, 0, 1, 0, 0]
    ])
    return np.concatenate((pat,)*2, axis=1)


def hollow_square(sz=3):
    if not sz % 2:
        sz += 1
    pat = npo((sz, sz))
    pat[1:-1, 1:-1] = 0
    return pat


def hs_1():
    sz = 11
    pat = hollow_square(sz)
    bar = npz((1, sz))
    #bar[:, (0, sz-1)] = 1
    pat = np.concatenate((bar, pat, bar), axis=0)
    return pat


def misc1():
    return npa([
        [0, 0, 1, 0, 0],
        [0, 1, 1, 1, 0],
        [1, 1, 0, 1, 1],
    ])


def pre_long():
    return npa([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0],
        [0, 1, 0],
        [1, 0, 1],
    ])


def four_pi():
    q = np.pad(pre_long(), 1)
    q_ud = np.flipud(q)
    zerosw, zerosh = map(lambda y: y-0, q.shape)
    zeros_sep = npz((zerosw, zerosh), dtype=bool)

    vert = np.concatenate((q_ud, zeros_sep, q), axis=0)
    print(vert.dtype)

    horz = np.rot90(vert)
    out_shape = vert.shape[0], horz.shape[1]
    print(q.shape, q_ud.shape, horz.shape, vert.shape,
          zeros_sep.shape, out_shape, sep='\n')
    out = npz(out_shape)
    vert_j = out_shape[0]-horz.shape[1]
    horz_i = out_shape[1]-vert.shape[0]
    # out[:, vert_j:vert_j+q_ud.shape[0]]
    return vert, horz


def long_pattern():
    pat = npz((12, 7))
    pat[(0, 10), 1:4] = 1
    pat[1:4, (0, 4)] = 1
    pat[7:10, (0, 4)] = 1
    pat[(11, ), 2] = 1
    return pat


def long_pattern2():
    pat = npz((12, 3))
    pat[0, (0, 2)] = 1
    pat[1, (0, 2)] = 1
    pat[2, (0, 2)] = 1
    pat[(3, 4), 1] = 1
    #pat[8:11, (0, 4)] = 1
    #pat[(11, 15), 2] = 1
    return pat


def dunno():
    return npa([
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
        [1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1],
        [0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
        [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0]
    ])


def small_test():
    return npa([
        [1, 1, 0, 1],
        [1, 1, 0, 1],
        [0, 0, 0, 1],
        [1, 1, 1, 0]
    ])


def init_universe_pulsar():
    pat = npz((15, 15))
    line = npz(15)
    line[3:6] = 1
    line[9:12] = 1
    for ind in (1, 6, 8, 13):
        pat[ind] = line
        pat[:, ind] = line
    return pat


def random_symmetric_quad(sz=256, density=0.35):
    """
    Square grid of ``sz`` pixels in size (== length of side) with four-fold
    symmetry, meaning that the four quadrants of the square are actually
    reflections of a single quadrant about the horizontal and vertical midlines
    of the square. The quadrant itself is randomly filled with live cells of
    the given density, which is a number from 0 to 1; a density of 0 results in
    an empty quadrant, whereas a density of 1 in a completely filled quadrant.

    :param sz: Length of side of the square, in pixels.
    :param density: Proportion of live cells in each quadrant, as a number
        from 0 to 1.
    :return: The pattern as a NumPy unsigned byte array.
    """
    pat, hsz = npz((sz, sz)), sz // 2 - 0
    pat[:hsz, :hsz] = (npr.random((hsz, hsz)) < density).astype(ui8)
    pat += np.flipud(np.fliplr(pat))
    pat += np.flipud(pat)
    return pat


def random_sparse_diag1(sz=256, density=0.7):
    hsz = sz//2
    sz = hsz*2
    pat = random_symmetric_quad(sz=sz, density=1)
    hrr = npr.sample(hsz) > density
    rr = npz(sz)
    rr[:hsz] = hrr
    rr += rr[::-1]
    dr = np.diag(rr)
    dr += np.fliplr(dr)
    return (np.logical_and(pat, np.logical_not(dr))).astype(ui8)


def random_sparse_diag3(sz=256, n=2, padding=1):
    hsz = (sz - 2*n*padding) // 2
    sz = hsz*2
    pat = np.pad(random_symmetric_quad(sz=sz//n, density=1), padding)
    hpat = np.concatenate((pat, )*n, axis=1)
    pat = np.concatenate((hpat, )*n, axis=0)
    return pat


rsq, rsd1, rsd3 = (
    random_symmetric_quad,
    random_sparse_diag1,
    random_sparse_diag3,
)

# Initial conditions that produced some interesting runs. "Profiles" are just
# groupings of keywords to be passed to Conway() -- see the documentation for
# Conway.default_options for details. Note, however, that the profiles are
# most likely out of date as they depended on implementation details that have
# changed by now. I did not have time to verify all of them, so they're
# provided "as-is"...
profiles = (
    dict(  # 0
        shape=(256, ) * 2,
        init_pattern=(rsq, 200, 0.5),
        seed=1911858326,
    ),
    dict(  # 1
        # very cool -- it evolves membrane receptors at the end! (at the steady
        # state, which is a cycle of T=2)
        shape=(512, ) * 2,
        init_pattern=(rsq, 216, 1),
    ),
    dict(  # 2
        # has that T=30 oscillator at the steady state
        shape=(512,) * 2,
        init_pattern=(rsq, 200, 1),
    ),
    dict(  # 3
        shape=(512, ) * 2,
        init_pattern=(rsq, 510, 1),
        seed=1078576960,
    ),
    dict(  # 4
        shape=(128*7, ) * 2,
        init_pattern=(rsd3, 512+17*8, 4),  # 512+2*8, 12-14*8 !!
        # seed=952878327,
    ),
    dict(  # 5
        shape=(128*7, ) * 2,
        init_pattern=(rsd3, 512 + 19 * 8, 4),  # 512+2*8, 12-14*8 !!
        # seed=952878327,
    ),
    dict(  # 6
        shape=(128*7, ) * 2,
        init_pattern=(rsd3, 512 + 48 * 8, 8),  # 512+2*8, 12-14*8 !!
        # seed=952878327,
    ),
    dict(  # 7
        # so far 19-22 are great!
        shape=(128*7, ) * 2,
        init_pattern=(rsd3, 512 + 58 * 8, 4),  # 512+2*8, 12-14*8 !!
        # seed=952878327,
    ),
    dict(  # 8
        # great, and 9600 gens till equilibrium!
        shape=(128*11, ) * 2,
        init_pattern=(rsd3, 1 << 10, 2),  # 512+2*8, 12-14*8 !!
        # seed=952878327,
    ),
    dict(  # 9
        # another T=30 and great!
        shape=(128*11, ) * 2,
        # good: 512/4/10
        # great: 512/4/16, 512/4/9, 512/3/2, 512/3/30 (!), 512/3/10,
        #  512/2/3, 512/2/4, 512/2/17, 512/5/6, 1024/1/60
        # long: 512/4/13 (5902 frames), 512/4/7
        init_pattern=(rsd3, 512*2, 1, 68),
        # seed=952878327,
    ),
    dict(   # 10
        # Maybe my fav so far
        shape=(512, ) * 2,
        init_pattern=(rsd3, 1460, 12, 10),  # 512+2*8, 12-14*8 !!
        # seed=952878327,
    ),
)
