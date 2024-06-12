import rectangle_packing_solver as rps

# Define a problem
rectangles=[
    {"width": 1, "height": 2, "rotatable": True},
] * 5
problem = rps.Problem(rectangles)
print("problem:", problem)

# [Other Usages]
# We can also give a solution width (and/or height) limit, as well as progress bar and random seed
print("\n=== Solving with width/height constraints ===")
solution = rps.Solver().solve(problem=problem, height_limit=3, width_limit=4, show_progress=True, seed=1111)
print("solution:", solution)
rps.Visualizer().visualize(solution=solution)
