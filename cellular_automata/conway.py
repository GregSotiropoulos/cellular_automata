# -*- coding: utf-8 -*-

"""
.. module:: conway
   :platform: Unix, Windows
   :synopsis: Conway's Game of Life

.. moduleauthor:: Greg Sotiropoulos <greg.sotiropoulos@gmail.com>

Cellular Automata (CA) application. Conway's Game of Life is the "prototypical"
CA and is supported by default, although other rule systems can be easily
incorporated (see the :class:`ConwayFilters` documentation). The module can be
used to compute the evolution of a grid from some initial configuration.

|
Grids are NumPy 2D byte arrays that can only have 0s and 1s, representing
dead and live cells, respectively. Note that grids are "closed", wrapping
around in both dimensions. Thus what is south of the bottom row is the top row
(and what is north of the top row is the bottom row), and same goes for
left/right. Essentially the CA world lives on a torus.

|
When run, the application presents a GUI that allows one to use mouse and
keyboard shortcuts to view the evolving grid, pause it, draw or delete new
cells, save the current grid as an image, pickle the current :class:`Conway`
instance to save the entire simulation and more. See the :class:`ConwayGui`
documentation for details. In the GUI, you can press 'h' to bring up a help
dialog that lists the available keyboard and mouse shortcuts.

|
To use the module in client code without the GUI, it is sufficient to import
:class:`Conway`.
"""

from operator import attrgetter
from pathlib import Path
from types import SimpleNamespace as Sns
from itertools import repeat, islice
from time import perf_counter as t, strftime
from collections.abc import Sequence
from numbers import Number
from tkinter import PhotoImage
from tkinter.messagebox import askyesno as tk_yesno, \
    showerror as tk_showerror, showinfo as tk_showinfo
from tkinter.filedialog import asksaveasfilename as tk_savefile, \
    askopenfilename as tk_openfile
import os
import pickle
# import logging

# non-stdlib modules
import numpy as np
from numpy import random as npr, uint8 as ui8, intp, ndarray as nda
from scipy.ndimage import convolve
import keyboard as kb
import matplotlib as mpl
# I know the linter doesn't like this but apparently use() has to be called
# before importing pyplot -- see https://stackoverflow.com/a/56320309/8441268
mpl.use('TkAgg')
from matplotlib import (
    pyplot as plt,
    backend_bases as mpl_bb,
    backend_tools as mpl_tools,
)
from matplotlib.animation import FuncAnimation as FuncAni

# package module
import cellular_automata.conway_patterns as patterns

# logger = logging.getLogger(__name__)
__docformat__ = 'reStructuredText'
__author__ = 'Greg Sotiropoulos <greg.sotiropoulos@gmail.com>'
__version__ = 1, 0, 1
__all__ = 'Conway', 'ConwayGui'


class ConwayFilters:
    """
    Filters (ie convolution matrices) used to compute the next generation. A
    ``ConwayFilters`` instance is basically a mapping (like a dictionary) that
    contains two built-in filters (which you cannot delete). User-defined
    filters can be added/deleted that implement different CA update rules.

    The filters contained in :class:`ConwayFilters` instances are 2-tuples. In
    such a tuple ``p``, ``p[0]`` is the 2D convolution array (typically 3x3)
    and ``p[1]`` is a sequence of numbers. At each time-step, the current grid
    is first filtered (convolved) with ``p[0]``. The grid is then updated as
    follows (in pseudo-Python):

        ``grid[x, y] = filtered_grid[x, y] in p[1]``

    for every (x, y) coordinate. This results in a new binary grid and
    the process repeats.
    """
    _built_in = 'conway', 'conway5'

    def __init__(self):
        c3, c5 = self._built_in
        fd = self._filters = {}

        # 'conway' filter, corresponding to the original update rule defined
        # by Conway in Game of Life
        filt = 2 + np.zeros((3, 3), dtype=ui8)
        filt[1, 1] = 1
        # if the value of a pixel after convolution is in the following
        # sequence, the pixel is set to 1, otherwise to 0.
        set_alive_vals = 6, 5, 7
        fd[c3] = filt, set_alive_vals

        # 'conway5' filter -- like the original but neighbourhood is +/-2
        # cells from current (on either axis), instead of +/-1.
        filt = np.zeros((5, 5), dtype=ui8)
        filt[2, 0] = filt[2, 4] = filt[0, 2] = filt[4, 2] = 1
        filt[1:-1, 1:-1] = 2
        filt[2, 2] = 1
        set_alive_vals = 8, 7, 5
        fd[c5] = filt, set_alive_vals

    def __contains__(self, filt_name):
        return filt_name in self._filters

    def __getitem__(self, filt_name):
        return self._filters.__getitem__(filt_name)

    def __setitem__(self, filt_name, filt):
        if filt_name in self._built_in:
            raise ValueError(f'cannot modify built-in filter "{filt_name}".')
        filt_array, truthy_range = filt
        if not (isinstance(filt_array, nda) and
                isinstance(truthy_range, Sequence) and
                all(map(isinstance, truthy_range, repeat(Number)))):
            raise TypeError(
                'filter must be a tuple of the form (filt_array, truthy_range)'
                f' with filt_array an instance of {nda.__name__} and '
                'truthy_range a sequence of numbers. Cell will be set alive '
                '(1) if the filtered version is one of those numbers, '
                'otherwise it will be set to dead (0).'
            )

    def __delitem__(self, filt_name):
        if filt_name in self._built_in:
            raise ValueError(f'cannot remove built-in filter "{filt_name}".')
        self._filters.pop(filt_name, None)


class Conway:
    """
    The core class of the module. It represents a cellular automaton (CA)
    simulation as a sequence of grids (two-dimensional NumPy uint8 arrays).

    ``Conway`` instances implement a read-only subset of both the Sequence and
    the Mapping interface. This means that grids in a simulation can be
    accessed both with a numeric or slice index and a 'grid' index.
        - In the case of a numeric index, it is straightforward: ``inst[i]``
          retrieves the grid after i generations, as a NumPy array. As in lists
          and tuples, ``i`` can be negative, with the same semantics.
        - In the case of a slice, a lazy (``map``-based) iterator is returned,
          eg ``inst[i:i+4]``. This is done to conserve memory and is thus
          unlike list slices, which produce new lists in a non-lazy way.
        - In the case of a grid, ``inst[grid]`` is the numeric index (ie the
          position) of the grid in the sequence computed so far. [Yes, NumPy
          arrays are not hashable and thus cannot be used as keys; this is
          made possible by transparently encoding/decoding the array into/from
          a ``bytes`` object.]

    The read-only constraint is to provide guarantees that it represents a
    valid history (ie that each frame comes from the immediately previous frame
    through the application of the CA rule defined in the current Filter).

    The following are all valid ways of indexing a ``Conway`` instance:

    >>> cw = Conway()
    >>> i, grid = cw.get_grid(5)
    >>> i == 5
    True

    >>> np.array_equal(cw[-2], cw[4])
    True

    >>> cw[cw[3]]
    3

    Please read the ``Conway.default_options()`` code for keywords that you can
    use in the constructor to set simulation options.
    """

    @staticmethod
    def is_pow2(x):
        """
        Determines if an integer (or an array of integers -- fully vectorized)
         is a power of 2.

        :param x:
        :return:
        """
        # x = np.array(x); return np.logical_and(x, np.logical_not(x & (x-1)))
        return 1 << np.log2(x).astype(intp) == x

    @staticmethod
    def round_to_pow2(x):
        """
        Round a number (or array of numbers) to closest power of 2.

        :param x: Number to be rounded.
        :return: Resulting integer power of 2.
        """
        return 1 << np.round(np.log2(x)).astype(intp)

    @staticmethod
    def is_mult8(x):
        """
        Determines if a number (or array of numbers) is a multiple of 8.

        :param x: Number to be checked.
        :return: Boolean result.
        """
        return not np.mod(x, 8).any()

    @staticmethod
    def round_to_mult8(x):
        """
        Round a number (or array of numbers) to closest multiple of 8.

        :param x: Number to be rounded.
        :return: Resulting integer multiple of 8.
        """
        return np.round(x / 8).astype(intp) << 3

    @classmethod
    def default_options(cls):
        """
        Class method that returns the default options of Conway instances.

        :return: A SimpleNamespace instance, in which each option can be
            accessed as an attribute.
        """
        return Sns(
            # name of filter to use for grid update (see ConwayFilters
            # documentation)
            filter_name='conway',

            # rng (random number generator) seed (for reproducible results
            # when a initial pattern is randomly generated. See module
            # conway_patterns for randomized patterns that you could use to
            # initialize your grid.
            seed=0,

            # H x W; keep both a power of 2 or a multiple of 8 -- see next
            shape=(512//2, ) * 2,

            # maximum number of cells -- you should not attempt grids bigger
            # than this in general as it will be very resource-consuming
            max_size=4096**2,

            # rounding method, when H and/or W are not multiple of 8/power of 2
            round_to=cls.round_to_mult8,

            # Initial conditions -- a tuple consisting of the function used
            # to generate the initial pattern, followed by any parameters
            # that this function accepts. The pattern is then centered on the
            # an empty grid of the specified shape. For example, if
            # ``options.init_pattern`` is set to ``(f, 1, 'two')`` then
            # ``f(1, 'two')`` would be called and its return value (a NumPy
            # array) would be placed on an empty grid (another Numpy array of
            # zeros of shape ``opts.shape``). Note that f() must return an
            # array that is no larger in either dimension than the grid.
            init_pattern=(patterns.random_symmetric_quad, 200, 1),

            # file name mask (for input in time.strftime) for save operations
            # in Conway (saving the entire instance data) or ConwayGui
            # (calling Conway.save() or saving the current image to a PNG file)
            save_fn_fmt='conway_%Y-%m-%d_%H.%M.%S',
        )

    @classmethod
    def _new_grid(cls, **opts):
        opts = Sns(**{**vars(cls.default_options()), **opts})
        f, *args = opts.init_pattern
        return cls._center(f(*args), np.zeros(opts.shape, dtype=ui8))

    @staticmethod
    def _center(pattern, grid):
        psh = pattern.shape
        dh, dw = (grid.shape - np.array(psh))//2
        if dh < 0 or dw < 0:
            raise ValueError('grid must be at least as big as the '
                             '(bounding box of) the pattern!')
        grid[dh:dh+psh[0], dw:dw+psh[1]] = pattern
        return grid

    @classmethod
    def _validate_shape(cls, x, options=None):
        try:
            x = np.asarray(x)
        except RuntimeError:
            raise TypeError(
                'shape must be a Sequence of exactly two integers or a NumPy '
                'numeric 1D array of exactly two elements'
        )
        if x.ndim != 1:
            raise ValueError(
                'shape must be a vector (tuple, list or 1D NumPy array)'
            )
        if x.size != 2:
            raise ValueError('shape must have exactly two elements')

        if options is None:
            options = cls.default_options()

        max_size = options.max_size
        if np.prod(x) > max_size:
            raise ValueError(
                f'grid size (number of cells) must be at most {max_size}'
            )
        return tuple(x)

    @classmethod
    def _validate_grid(cls, x):
        try:
            x = np.asarray(x).clip(0, 1).astype(ui8)
        except RuntimeError:
            raise TypeError(
                'input must be a NumPy numeric 2D array or a sequence that '
                'can be converted to one'
            )
        return x

    def __init__(self, *args, **kwa):
        cls = __class__
        self._options = opts = cls.default_options()
        opts_dic = vars(opts)
        opts_dic.update(kwa)
        if opts.seed is not None:
            npr.seed(opts.seed)

        # if there is at least one (non-keyword) argument, it must be either a
        # NumPy array or another Conway instance.
        if args:
            a = args[0]
            init_grid = (
                a[0] if isinstance(a, cls) else
                self._validate_grid(a)
            )
            opts.shape = init_grid.shape
        else:
            init_grid = self._new_grid(**opts_dic)
            # check that the grid shape conforms to the limits set in options
            opts.shape = self._validate_shape(init_grid.shape, opts)

        # this is stored to slightly speed up self.unpack()
        self._shape = opts.shape

        # filter to be used for state update -- see ConwayFilters documentation
        opts.filter = ConwayFilters()[opts.filter_name]

        # indicator that a cycle (ie a repeating sequence of frames) has been
        # detected. If not None (meaning no cycle yet), then it stores the
        # index of the first grid in the cycle. Used by self._evolve/get_grid
        self._cycle = None
        self._dic = {}
        self._add(init_grid)
        # Conway instances can be indexed by NumPy arrays (by value), integers
        # slices and tuples. The latter are interpreted as arguments to slice()
        self._getitem = {
            nda: self._gi_nda,
            int: self._gi_int,
            slice: self._gi_slice,
            tuple: self._gi_tuple
        }

    def __contains__(self, k):
        return isinstance(k, nda) and __class__.pack(k) in self.dic

    def __getitem__(self, k):
        return self._getitem[type(k)](k)

    def __len__(self):
        return len(self.dic)

    def __iter__(self):
        return map(self.unpack, self.dic)

    __setitem__ = __delitem__ = None

    def _gi_nda(self, grid):
        return self.dic[__class__.pack(grid)]

    def _gi_int(self, i):
        err = IndexError(f'grid sequence index out of range')
        try:
            if i < 0:
                i += len(self.dic)
                if i < 0:
                    raise err
            return next(self._gi_slice(slice(i, i+1)))
        except StopIteration:
            raise err

    def _gi_slice(self, s):
        return map(self.unpack, islice(self.dic, *s.indices(len(self))))

    def _gi_tuple(self, t):
        return self._gi_slice(slice(*t))

    def _add(self, grid):
        self.dic[__class__.pack(grid)] = len(self)

    @staticmethod
    def pack(grid):
        """
        Pack a grid and convert it to ``bytes``. This reduces storage
        requirements for the grid sequence and provides a hashable
        representation of the grid (a NumPy array that is not hashable).
        Packed arrays can thus be used "by value" as dictionary keys.

        :param grid: A 2D NumPy array representing the CA grid.
        :return: The packed grid.
        """
        return np.packbits(grid).tobytes()

    def unpack(self, packed_grid):
        """
        The inverse operation of Conway.pack(). Unlike the latter, this is an
        instance method as it depends on state (the shape of the unpacked grid)

        :param packed_grid: A ``bytes`` representation of a packed grid.
        :return: The grid as it is normally -- a 2D NumPy uint8 array.
        """
        h, w = self._shape
        packed_grid = np.frombuffer(packed_grid, dtype=ui8)
        return np.unpackbits(packed_grid)[:h*w].reshape(h, w)

    @property
    def dic(self):
        """
        A dictionary of computed grids, serving as a cache (and a memo, for
        detecting a possible cycle, ie a grid sequence that repeats forever).

        Keys in this dictionary are the ``bytes`` representations of (packed)
        grids. Values are just ordinals, starting from zero (you can think of
        them as number of generations). This all works because as of Python
        3.6, ``dict``s are ordered.

        :return: The dictionary instance.
        """
        return self._dic

    @property
    def options(self):
        return self._options

    def _evolve(self, gens=1):
        """
        Evolve the grid, ie compute successive generations. Computed grids
        are stored in an (ordered) dictionary, mainly so that cycles can be
        detected.

        :param gens: The number of generations. Grid is iteratively updated
            ``gens`` times.
        :return: The final grid, after ``gens`` generations from the current.
        """
        if not (np.isinf(gens) or isinstance(gens, int) and gens > 0):
            print('gens =', gens)
            raise TypeError('gens must be a positive integer or infinity.')

        opts = self.options
        filt, live_vals = opts.filter
        grid, n = self[-1], len(self)

        for k in range(n, n+gens):
            conv_grid = convolve(grid, filt, mode='wrap')
            # alternative convolution method from scipy.signal
            # conv_grid = convolve2d(grid, filt, mode='same', boundary='wrap')

            # compare the convolved grid against the "set-alive" list of values
            grid = np.isin(conv_grid, live_vals).astype(ui8)

            # if we've seen this grid before, it's a cycle
            if grid in self:
                if self._cycle is None:  # if this is the first recurring grid
                    self._cycle = c = self[grid]
                    print(f'Cycle of length {k-c} detected '
                          f'(frame {c} == frame {k})')
                break
            else:
                self._add(grid)
        return grid

    def get_grid(self, i):
        """
        Get grid in the ``i``th generation (ie after ``i-1`` updates).

        :param i: The generation number, eg ``i = 3`` would compute (or
            retrieve from cache) the first four generations.
        :return: A 2-tuple. The first element is the absolute (ie positive)
            index of the ``i``th grid. This will be equal to ``i`` if there is
            no cycle encountered yet, otherwise it may be less than ``i``. This
            can happen if the latter is greater than the cycle's last index.
            The second element is the ``i``-th grid as a NumPy array.
        """
        n = len(self)
        if i < n:
            return i, self[i]
        c = self._cycle
        if c:
            i = c + (i - c) % (n - c)
            return i, self[i]
        # i-th grid not in cache; compute all grids from current last to i-th
        return i, self._evolve(i - n + 1)

    def tofile(self, file_name=None):
        """
        Save the simulation, that is, the Conway instance that contains the
        sequence of grids computed so far.

        :param file_name: Name of file.
        """
        file_name = (
            strftime(self.options.save_fn_fmt) + '.conway'
            if file_name is None else
            file_name
        )
        with open(file_name, 'wb') as file:
            pickle.dump(self, file)

    @staticmethod
    def fromfile(file_name):
        """
        Load a CA run (== simulation) from a pickle file (one previously saved
        with ``inst.save()``).

        :param file_name: Name of pickle file.
        """
        with open(file_name, 'rb') as file:
            return pickle.load(file)


class ConwayGuiEvents:
    """
    Wrapper for mouse and keybord events in the GUI.
    """
    def __init__(self, gui):
        self._gui = gui
        # check if we're starting the animation in a paused state
        self.paused = gui.options.animation.start_paused
        self.paused_grid = self.imgrid.copy() if self.paused else None

        # mouse events setup

        # set up a namespace where m.R returns the code for the Right mouse
        # button, m.L for Left etc
        try:
            btns = mpl_bb.MouseButton.__members__
            self.mouse = m = Sns(**{name[0]: int(btns[name]) for name in btns})
        except AttributeError:
            # if there is an older matplotlib version that doesn't have the
            # MouseButton attribute (which is a simple IntEnum)
            self.mouse = m = Sns(L=1, M=2, R=3, B=8, F=9)

        # tracks which buttons are pressed
        m.down = np.zeros(len(vars(m)), dtype=ui8)

        # placeholder for selections using the mouse
        # (eg to copy part of the grid)
        m.selection = np.zeros((2, 2), dtype=intp) - 1

        # mouse and keyboard event callbacks
        cb = dict(
            button_press=self._on_mouse_down,
            button_release=self._on_mouse_up,
            scroll=self._on_wheel,
            motion_notify=self._on_mouse_move,
            key_press=self._on_key_press,
            close=self._on_close,
        )

        # event connections to the canvas
        connect = self.canvas.mpl_connect
        self._connections = [connect(k + '_event', v) for k, v in cb.items()]

    @staticmethod
    def info():
        return """
    Keyboard shortcuts
    ---------------------------------------------------------------------------
    h:			show this help dialog
    Space:			pause/resume animation
    Del:			clear grid *
    Left:			go back one generation *
    Right:			advance one generation *
    l (L):			load Conway instance from 
     			saved pickle file
    i:			restart animation from initial grid
    o:			load PNG image file as grid 
    s:			save grid to PNG (auto-generated 
    			filename)
    Ctrl+s:			save grid to PNG (filename 
    			selected via dialog)
    
    Mouse shortcuts
    ---------------------------------------------------------------------------
    Wheel up:		increase fps (frames per second)
    Wheel down:		decrease fps
    Left click:		erase pixel *, **
    Right click:		paint pixel *, **
    Shift + mouse move:	select and copy rectangle ***
    
    *	only when the animation is paused
    **	you can keep the left/right mouse button pressed 
    	while moving the mouse to erase/paint en masse
    ***	can be accessed via gui_inst.selection()
    """

    @property
    def gui(self):
        return self._gui

    @property
    def figure(self):
        return self.gui.fig

    @property
    def window(self):
        return self.gui.window

    @property
    def canvas(self):
        return self.figure.canvas

    @property
    def ev_src(self):
        return self.gui.ani.event_source

    @property
    def image(self):
        return self.gui.image

    @property
    def imgrid(self):
        """
        Retrieve grid from current figure.

        :return: The CA grid, in the form of a NumPy array.
        """
        return self.image.get_array().filled(fill_value=0)

    def _cleanup(self):
        cb_reg = mpl_tools.cbook.CallbackRegistry()
        # loop over truthy values only
        for cb in filter(None, cb_reg.callbacks):
            cb_reg.disconnect(cb)
        delattr(self, '_gui')

    def _mouse_event2ji(self, e):
        """
        Get coordinates (in pixels) under the mouse cursor from a mouse
        event.

        :param e: The mouse event, an instance of mpl.backend_bases.MouseEvent
        :return: A 2-element NumPy vector with the coordinates.
        """
        yx = e.ydata, e.xdata
        return (
            np.full(2, -1, dtype=intp)
            if None in yx else
            np.minimum(np.floor(yx), self.imgrid.shape).astype(intp)
        )

    # event callbacks

    def _on_mouse_down(self, e):
        m, btn = self.mouse, e.button
        m.down[:] = btn == m.down
        if kb.is_pressed('left shift'):
            if btn == m.L:
                ji = self._mouse_event2ji(e)
                if (ji >= 0).all():
                    m.selection[0, :] = ji
        self._on_mouse_move(e)

    def _on_mouse_up(self, e):
        m, btn = self.mouse, e.button
        m.down[:], sel = (btn == m.down), m.selection
        if btn == m.L:
            if kb.is_pressed('left shift'):
                sel[1, :] = self._mouse_event2ji(e) + 1
                if (sel >= 0).all():
                    grid = self.imgrid
                    x_min, y_min = np.min(sel, axis=0)
                    x_max, y_max = np.max(sel, axis=0)
                    sub_array = grid[x_min:x_max, y_min:y_max]
                    self.gui._selection = sub_array
                    print(repr(sub_array))
        sel[:] = np.full((2, 2), -1, dtype=np.int8)

    def _on_wheel(self, e):
        if not self.paused:
            gui = self.gui
            ani_opts = gui.options.animation
            wheel_delta = (
                (1 - 2*(e.button == 'up')) * ani_opts.wheel_step/100
            )
            ev_src = self.ev_src
            interval_prev = ev_src.interval
            interval = round(gui._clip_i(interval_prev * (1 + wheel_delta)))
            # this is necessary due to round-off issues: when interval is small
            # a 10% change may lead to the same integer after rounding/flooring
            if interval == interval_prev:
                interval += np.sign(wheel_delta)
            ev_src.interval = interval

    def _on_mouse_move(self, e):
        if not kb.is_pressed('left shift'):
            m, btn = self.mouse, e.button
            i, j = self._mouse_event2ji(e).flatten()
            lr = m.L, m.R
            if btn in lr and self.paused:
                # moving the mouse while paused:
                #    - kills cells if the left button is pressed
                #    - creates cells if the right button is pressed
                px_val = lr.index(btn)
                grid = self.imgrid
                grid[i, j] = px_val
                self.image.set_data(grid)
                self.canvas.draw_idle()

    def _on_key_press(self, e):
        gui, grid, key = self.gui, self.imgrid, e.key.lower()

        if key == 'h':
            tk_showinfo(
                parent=self.window,
                title='Game of Life help',
                message=__class__.info()
            )

        elif key == ' ':
            self.toggle_pause()

        # clear grid (kill all cells)
        elif key == 'delete':
            if self.paused:
                self.image.set_data(np.zeros_like(grid))

        elif key == 'l':
            w = self.window
            open_fn = tk_openfile(
                parent=w,
                initialdir=os.getcwd(),
                title='Load pickled Conway instance',
                filetypes=(('Conway pickle files', '*.conway'), )
            )
            if open_fn:
                self._load_new_grid(Conway.fromfile(open_fn))

        # go back to the initial frame
        elif key == 'i':
            grid, gui._frame = gui.conway[0], 0
            self.image.set_data(grid)
            if self.paused:
                self.paused_grid = grid.copy()

        # open image file
        elif key == 'o':
            w = self.window
            open_fn = tk_openfile(
                parent=w,
                initialdir=os.getcwd(),
                title='Open image file as initial grid',
                filetypes=(
                    ('PNG files', '*.png'),
                    ('All files', '*.*')
                )
            )
            if open_fn:
                new_grid = plt.imread(open_fn)
                # if image is RGB, average over the 3 channels (ignoring alpha)
                new_grid = (
                    # anything whiter than mid-gray is considered a live cell
                    np.mean(new_grid[:, :, :3], axis=2) > 0.5
                ).astype(ui8)
                try:
                    self._load_new_grid(Conway(new_grid))
                except RuntimeError(e):
                    tk_showerror(
                        parent=w,
                        title='Error loading image',
                        message=e
                    )

        # save image to file
        elif key in ('s', 'ctrl+s'):
            save_fn = strftime(gui.conway.options.save_fn_fmt) + '.png'
            if key != 's':  # if Ctrl+s is pressed
                save_fn = tk_savefile(
                    parent=self.window,
                    initialdir=os.getcwd(),
                    title='Save grid to image file',
                    filetypes=(('PNG files', '*.png'), )
                )
            if save_fn:
                if not save_fn.lower().endswith('.png'):
                    save_fn += '.png'
                plt.imsave(save_fn, grid, cmap=self.image.cmap)

        # move back/forward a generation
        elif key in ('left', 'right'):
            if self.paused:
                gui = self.gui
                sense = 2*len(key) - 9  # left -> -1, right -> 1
                frame = gui._frame + sense
                if frame >= 0:
                    gui._frame = frame
                    k, grid = gui.conway.get_grid(frame)
                    self.paused_grid = grid.copy()
                    self.image.set_data(grid)

        self.canvas.draw_idle()

    def _load_new_grid(self, conway):
        img = self.image
        img.set_data(conway[0])
        h, w = conway[0].shape
        # set extent again as image may have different dimensions
        img.set_extent((0, w, h, 0))
        if not self.paused:
            self.gui._new_animation(conway)

    def toggle_pause(self):
        """
        Hopefully self-explanatory; pauses the animation, or resumes it if it
        was already paused.
        """
        if self.paused:
            self.paused = False
            self._on_change()
            self.ev_src.start()
            self.paused_grid = None
        else:
            self.ev_src.stop()  # stop Animation
            self.paused = True
            self.paused_grid = self.imgrid.copy()

    def _on_change(self):
        imgrid, gui = self.imgrid, self.gui
        try:
            if (
                not np.array_equal(self.paused_grid, imgrid) or
                imgrid not in gui.conway
            ):
                # last displayed frame has been modified, so we spawn a new
                # "session"
                gui._new_animation(Conway(imgrid))
        except KeyError:
            gui._cleanup()
            raise RuntimeError('internal state is irrecoverably corrupt! '
                               'Please restart the application.')

    def _on_close(self, e):
        gui = self.gui
        save = gui.options.save_on_exit
        if save in (None, 'ask'):
            self.ev_src.stop()
            window = gui.window
            window.withdraw()
            if tk_yesno(
                title='Save simulation',
                message='Do you want to save the current simulation (the '
                        'Conway instance that contains the sequence of grids '
                        'computed so far) to a pickle file?',
                parent=window
            ):
                gui.conway.tofile()
        elif save and save != 'no':
            gui.conway.tofile()


class ConwayGui:
    """
    Graphical User Interface (GUI) for the Game of Life. Creating an instance
    of ``ConwayGui`` opens a figure window and sets up an animation of the
    grid's evolution.

    Please read the ``ConwayGui.default_options()`` code for keywords that you
    can use in the constructor to set GUI-related options. Pressing the 'h' key
    will display information on available mouse and keyboard commands -- make
    sure to have a look. Alternatively, the shortcut information can be printed
    by ``print(gui.events.info())`` where ``gui`` is a ``ConwayGui`` instance.

    In order for the user not to lose the current simulation data after making
    edits to the frame when paused, the application automatically saves the
    entire simulation (the sequence of grids that have been procuded thus far)
    in ``pickle`` format upon resuming the animation. Simulations can be
    "replayed" by loading a pickled ``Conway`` instance through the GUI.

    Apart from pickled instances, images can also be loaded, in which case the
    image is converted to grayscale, thresholded at mid-gray (meaning that any
    pixel with value > 0.5 is considered a live cell, otherwise a dead one) and
    the resulting 2D NumPy array is passed to the ``Conway`` constructor.
    """
    @staticmethod
    def _clip_i(interval):
        # clip animation interval between 1 ms and 1 min. This is mainly to
        # prevent zeros (and therefore divisions by zero) in fps/interval.
        return np.clip(interval, 5, 60_000)

    @staticmethod
    def _clip_f(fps):
        # the limits here correspond to the interval limits defined above
        return np.clip(fps, 1/60, 200)

    @classmethod
    def default_options(cls):
        """
        Default GUI options.

        :return: A SimpleNamespace instance, in which each option can be
            accessed as an attribute.
        """
        return Sns(
            # whether to save the simulation (the Conway instance) when exiting
            # the GUI. Possible values: 'yes'/1/True; 'no'/0/False; 'ask'/None
            save_on_exit='ask',

            # Whether to save the current simulation after editing a paused
            # frame with the mouse. Resuming the animation spawns a new
            # simulation with the edited frame as initial conditions.
            save_on_change=True,

            # animation-related options
            animation=Sns(
                start_paused=False,  # start the animation in a paused state
                fps=30,  # target frame rate
                wheel_step=10,  # fps percentage change on wheel scroll
            ),
        )

    def __init__(self, **kwa):
        """
        Constructor (initializer, to be precise) that only accepts keyword
        arguments.

        :param kwa: Keyword arguments of the same name as those in
            ``ConwayGui.default_options()`` are used to customize some aspects
            of the GUI. The rest of the keywords are passed to ``Conway()``
        """
        self._options = opts = self.default_options()
        ani_opts = opts.animation
        ani_opts.start_paused = kwa.pop('start_paused', ani_opts.start_paused)
        opts_dic = vars(opts)

        # extract any kwa that match (top-level) option names,
        # pass the rest of the keywords to the Conway constructor
        opts_dic.update((k, kwa.pop(k)) for k in opts_dic if k in kwa)
        conway = Conway(**kwa)

        self._figure = plt.figure(f'Game of Life', frameon=False)
        self._image = img = plt.imshow(
            conway[0], cmap='gray', aspect='equal', animated=True
        )
        h, w = conway._shape
        img.set_extent((0, w, h, 0))
        for cmd in 'grid grid_minor save xscale yscale zoom'.split():
            plt.rcParams['keymap.' + cmd].clear()
        plt.axis('off')
        plt.tight_layout(pad=0)
        self._window = window = plt.get_current_fig_manager().window

        # application window icon
        iconfile = Path(Path(__file__).parent, 'glider.png')
        if iconfile.exists():
            icon = PhotoImage(file=str(iconfile), master=window)
        window.tk.call('wm', 'iconphoto', window._w, icon)

        self._events = ConwayGuiEvents(self)

        # this will store the selection, in the form of a NumPy array,
        # made by dragging the mouse while holding down the Shift key
        self._selection = None
        self._new_animation(conway)

        # This is necessary when running the app outside of IPython. Even
        # interactive mode needs it. It also needs to be at the end (after
        # all other matplotlib code has been called)
        plt.show()

        # quick and dirty figure height calculation -- it could be quite off
        # depending on your desktop configuration, number of monitors and OS
        scr_h, scr_w = window.winfo_screenheight(), window.winfo_screenwidth()
        fig_h = int(scr_h * 0.9)
        fig_w = int(min(fig_h * w / h, scr_w * 0.9))
        fig_x, fig_y = (scr_w - fig_w)//2, 0
        window.geometry(f'{fig_w}x{fig_h}+{fig_x}+{fig_y}')

    def _new_animation(self, conway):
        """
        Starts a new animation.
        """
        opts = self.options
        if opts.save_on_change and hasattr(self, '_conway'):
            self.conway.tofile()

        self._conway = conway

        # used for tracking the frame currently being displayed, for the case
        # where left/right arrows have been used to move
        self._frame = 0

        if hasattr(self, '_animation'):
            ani = self.ani
            interval = ani.event_source.interval
            ani._stop()
        else:
            interval = 1000/self._clip_f(opts.animation.fps)

        def frames_gen():
            """
            Generator for FuncAnimation's ``frames`` argument. The next 'frame'
            is actually the current timestamp. There is a frame counter (or,
            more precisely, index) kept in ``self._frame`` and incremented with
            every call to the animations update function (``self._update``).
            """
            while True:
                yield t()

        self._animation = a = FuncAni(
            self.fig,
            self._update,
            frames=frames_gen,
            interval=interval,
            # modifying the following options will generally slow things down
            repeat=False,
            blit=False,
            cache_frame_data=False,
        )
        self._t0 = self._t1 = t()
        print(f'Started animation (seed={conway.options.seed})')

    def _update(self, t1):
        """
        Animation update callback (second argument in FuncAnimation). It is
        called every time the animation's timer fires (see mpl's TimedAnimation
        documentation).

        :param t1: The current timestamp (as returned by time.perf_counter)
        :return: A one-element tuple containing the image of the new grid.
        """
        try:
            conway, events, ani = attrgetter('conway', 'events', 'ani')(self)
        except AttributeError:
            return self._cleanup()
        if events.paused:
            ani.event_source.stop()
            return ()
        im = self.image
        self._frame = frame = self._frame + 1
        i, grid = conway.get_grid(frame)
        t2 = t()
        im.set_data(grid)
        fps_mean = frame/(t2-self._t0)  # average fps
        dt, self._t1 = t1-self._t1, t1
        fps_inst = 1/dt  # (instantaneous) fps
        print(
            f'Frame {i} computed in {1000*(t2-t1):.1f} ms, '
            f'drawn in {1000*dt:.1f} ms '
            f'(fps: {fps_inst:.1f}, mean fps: {fps_mean:.1f})',
            end='\r'
        )
        return im,

    @property
    def options(self):
        """
        Get GUI options.

        :return: The options, as a SimpleNamespace instant.
        """
        return self._options

    @property
    def conway(self):
        """
        Access the Conway instant associated with the GUI window.

        :return: Instant of Conway class.
        """
        return self._conway

    @property
    def window(self):
        return self._window

    @property
    def fig(self):
        return self._figure

    @property
    def image(self):
        return self._image

    @property
    def ani(self):
        return self._animation

    @property
    def events(self):
        return self._events

    @property
    def selection(self):
        """
        Get rectangular subarray from figure by left-clicking, dragging and
        releasing the mouse button, all while holding down Shift.

        :return: A NumPy array corresponding to the rectangle drawn.
        """
        return self._selection

    def _cleanup(self):
        print('Starting cleanup...', end=' ')
        for attr, meth in (
            ('_events', ConwayGuiEvents._cleanup),
            ('_animation', FuncAni._stop)
        ):
            try:
                delattr(self, attr)
                meth()
            except AttributeError:
                pass
        print('Cleanup complete.')
        return ()


def main():
    """
    Called when the module is run. Small wrapper for ``ConwayGui()``.

    :return: A ConwayGui instance.
    """
    pidx = 1  # pick profile #1 from conway_patterns.profiles
    profile = patterns.profiles[pidx]
    # The GUI will open in a paused state. Hit the space key to start the
    # animation. Alternatively, start_paused can be set to False (the default).
    return ConwayGui(**profile, start_paused=True)


if __name__ == '__main__':
    # conway_gui contains all information a user might need, such as the
    # associated Conway instance
    conway_gui = main()
