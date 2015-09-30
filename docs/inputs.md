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

# Serial output for sample inputs
Running times were taken on a Phase 3 node.
- input_128x128.params
```
./lbm -a ./av_vels.dat -f ./final_state.dat -p ../inputs/input_128x128.params
==done==
Reynolds number:		2.039107350801E+01
Elapsed time:			103.457805 (s)
Elapsed user CPU time:		103.205310 (s)
Elapsed system CPU time:	0.004999 (s)
```

- input_pipe.params
```
./lbm -a ./av_vels.dat -f ./final_state.dat -p ../inputs/input_pipe.params
==done==
Reynolds number:		2.889654348141E+01
Elapsed time:			13.158120 (s)
Elapsed user CPU time:		13.123005 (s)
Elapsed system CPU time:	0.002999 (s)
```

- input_rect.params
```
./lbm -a ./av_vels.dat -f ./final_state.dat -p ../inputs/input_rect.params
==done==
Reynolds number:		3.740889367644E+01
Elapsed time:			212.361721 (s)
Elapsed user CPU time:		211.864791 (s)
Elapsed system CPU time:	0.001999 (s)
```
