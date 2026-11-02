import matplotlib.pyplot as plt
import numpy as np

# Create a figure
fig, ax = plt.subplots(figsize=(8, 6))

# Draw a line segment
x = [0, 3, 5, 2]
y = [0, 2, 1, 3]
ax.plot(x, y, 'gray', alpha=0.5, linewidth=2, label='Original line')

# Add arrows on each segment of the line
for i in range(len(x)-1):
    # Calculate the midpoint of the segment
    mid_x = (x[i] + x[i+1]) / 2
    mid_y = (y[i] + y[i+1]) / 2
    
    # Calculate the direction vector
    dx = x[i+1] - x[i]
    dy = y[i+1] - y[i]
    
    # Add arrow at the midpoint
    ax.annotate('', xy=(mid_x + dx/4, mid_y + dy/4), 
                xytext=(mid_x - dx/4, mid_y - dy/4),
                arrowprops=dict(arrowstyle='->', color=f'C{i}', 
                              lw=2, mutation_scale=15))

# Add markers at the endpoints of the line segment
ax.scatter(x, y, c='red', s=100, zorder=5, label='Endpoints')

# Set figure properties
ax.set_xlim(-1, 6)
ax.set_ylim(-1, 4)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.legend()
plt.title('Arrow Example on Line Segment')
plt.show()