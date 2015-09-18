# Input file

This is an example input file, with comments for documentation. *DO NOT* include the comments in your real input file.

```
128                # number of cells in x direction
128                # number of cells in y direction
40000              # number of timesteps
10                 # Reynolds dimension
0.1                # fluid density
0.005              # acceleration
1.85               # omega
accelerate row 4   # instruction to accelerate row/column i
6 obstacles        # n obstacles
0 0 100 1          # list of obstacle positions, see below
0 0 1 100
0 99 100 100
99 0 100 100
80 20 81 80
10 40 40 41

```

The obstacles are defined with two co-ordinates defining the bottom-left and top-right of the obstacle. The co-ordinates range from 0 to 100. This is independent of grid size. For example, `0 0 100 1` defines an obstacle along the bottom of the grid, from the bottom left (0,0) to the bottom right (100,1).
