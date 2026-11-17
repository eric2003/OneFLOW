import matplotlib.pyplot as plt
import numpy as np

# Example 4: How to truly center a line
plt.figure(figsize=(6, 6))

# Method 1: Manually set symmetric coordinate ranges
plt.plot([-1, 1], [0, 0], 'r-', linewidth=3, label='Horizontal line y=0')
plt.plot([0, 0], [-1, 1], 'b-', linewidth=3, label='Vertical line x=0')

# Key: Set symmetric axis limits
plt.xlim(-2, 2)
plt.ylim(-2, 2)

plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)  # x-axis
plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)  # y-axis
plt.grid(True)
plt.title('Centered by setting symmetric ranges')
plt.legend()
plt.gca().set_aspect('equal')  # Equal aspect ratio
plt.show()