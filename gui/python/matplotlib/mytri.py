import matplotlib.pyplot as plt
import numpy as np

xy = np.asarray([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.5, 0.866]
    ])

x = xy[:, 0]
y = xy[:, 1]

triangles = np.asarray([
    [0, 1, 2], ])
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.triplot(x, y, triangles )
ax.set_title('OneFLOW single triangle Test')

plt.show()
