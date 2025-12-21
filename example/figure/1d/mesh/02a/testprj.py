import matplotlib.pyplot as plt
import numpy as np

# Example 5: Canvas coordinates vs data coordinates
fig = plt.figure(figsize=(8, 6))

# Center in canvas coordinates (not data coordinates!)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # Graphics area centered in canvas

# Now draw lines in the graphics area
x_center = 5
y_center = 5
length = 2

# Horizontal line
plt.plot([x_center - length/2, x_center + length/2], 
         [y_center, y_center], 'r-', linewidth=3, label='Horizontal line')

# Vertical line
plt.plot([x_center, x_center], 
         [y_center - length/2, y_center + length/2], 'b-', linewidth=3, label='Vertical line')

# Set symmetric ranges to center line in graphics area
plt.xlim(x_center - length, x_center + length)
plt.ylim(y_center - length, y_center + length)

plt.grid(True)
plt.axhline(y=y_center, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=x_center, color='gray', linestyle='--', alpha=0.5)
plt.title('Line centered in graphics area')
plt.legend()
plt.show()