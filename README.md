# First approach

1. Add to each needed cut the half of the cut width.

2. Use this repo https://github.com/kotarot/rectangle-packing-solver
3. Try to find different solutions and use the waste as the chi2 to optimize the
   problem
4. Have a list with all the needs.
5. Shuffle the list and try sequentially put them in a board. When no more can
   be placed, start a new board. Continue until all the needs are placed.
6. Compute all the wastes. If it is the minimum found, store the configuration.
7. Continue with more shuffles (ideally all the configurations could be tested)

