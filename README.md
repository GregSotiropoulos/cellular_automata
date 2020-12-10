# Cellular Automata

Cellular Automata (CA) application. The core class, where the simulation is
performed is ``Conway``, named after John Conway's (RIP) Game of Life. This 
is the "prototypical" CA and is supported by default, although other rule 
systems can be easily incorporated (see the ``ConwayFilters`` documentation). 
The module can be used to compute the evolution of a grid from some initial 
configuration.

When run, the application presents a GUI that allows one to use mouse and
keyboard shortcuts to view the evolving grid, pause it, draw or delete new
cells, save the current grid as an image, pickle the current ``Conway`` 
instance to save the entire simulation and more. See the ``ConwayGui`` 
documentation for details. In the GUI, you can press 'h' to bring up a help 
dialog that lists the available keyboard and mouse shortcuts.

To use the module in client code without the GUI, it is sufficient to import
the ``Conway`` class.