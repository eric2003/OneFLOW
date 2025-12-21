import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(6, 6))

x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6], dtype=np.float64)

# Method 1: Manually set symmetric coordinate ranges
plt.plot([-1, 1], [0, 0], 'r-', linewidth=3, label='Horizontal line y=0')
plt.plot([0, 0], [-1, 1], 'b-', linewidth=3, label='Vertical line x=0')

# Key: Set symmetric axis limits
plt.xlim(-6, 6)
plt.ylim(-2, 2)

plt.title('Centered by setting symmetric ranges')
plt.legend()
plt.axis('equal')
#plt.axis('off')
plt.show()