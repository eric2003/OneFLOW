# Two-block laminar flat plate with ten geometric line segments

This case keeps the two-block structured topology while splitting each of the three vertical topological edges at y=0.02. The six vertical segments and four horizontal segments produce ten geometric lines. Each split pair maps back to one topological side.

The first 40 intervals between y=0 and y=0.02 are uniform, giving 41 points on each lower vertical segment. The upper vertical segments contain 45 points; together with the shared layer point this preserves 85 points on each complete vertical side.

Run this case through the normal OneFLOW grid task. The case is configured to generate a Plot3D grid and boundary file under `results/`, then convert the grid to OneFLOW format and write a grid-info summary.
