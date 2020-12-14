# Cellular Automata

[Cellular Automata](https://en.wikipedia.org/wiki/Cellular_automaton) (CA) application. 
The core class, where the simulation is performed, is `Conway`, named after
[John Horton Conway's](https://en.wikipedia.org/wiki/John_Horton_Conway) (RIP) 
[Game of Life](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life). This 
is the "prototypical" CA and is supported by default, although other rule 
systems can be easily incorporated (see the `ConwayFilters` documentation). 
The module can be used to compute the evolution of a grid from some initial 
configuration.

When run, the application presents a GUI that allows one to use mouse and
keyboard shortcuts to view the evolving grid, pause it, draw or delete new
cells, save the current grid as an image, pickle the current `Conway` 
instance to save the entire simulation and more. See the `ConwayGui`
documentation for details. In the GUI, you can press `h` to bring up a help 
dialog that lists the available keyboard and mouse shortcuts.

## Installation

In a terminal/command prompt, type:
```shell
python -m pip install git+https://github.com/GregSotiropoulos/cellular_automata.git
```

## Usage

To run the GUI application, type:
```commandline
python -m cellular_automata.conway
```
\
To use the module in client code without the GUI, it is sufficient to import `Conway`: 
```python
from cellular_automata.conway import Conway
```