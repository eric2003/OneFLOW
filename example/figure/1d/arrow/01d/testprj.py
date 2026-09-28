import matplotlib.pyplot as plt
import numpy as np

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Define line segment coordinates
x_start, y_start = 0, 0      # Starting point
x_end, y_end = 5, 3          # Ending point

# Draw the original line segment
ax.plot([x_start, x_end], [y_start, y_end], 
        'gray', alpha=0.5, linewidth=3, label='Original line segment')

# Calculate direction vector
dx = x_end - x_start
dy = y_end - y_start

# Calculate normalized direction vector
length = np.sqrt(dx**2 + dy**2)
dx_norm = dx / length
dy_norm = dy / length

# Define positions where to place arrows along the line
# Using fractions of the line length
arrow_positions = [0.2, 0.5, 0.8]  # At 20%, 50%, and 80% of the line

# Add arrows at specified positions
for i, fraction in enumerate(arrow_positions):
    # Calculate arrow position
    arrow_x = x_start + fraction * dx
    arrow_y = y_start + fraction * dy
    
    # Calculate arrow endpoints (start and end of arrow)
    arrow_length = 0.8  # Length of arrow relative to line
    
    # Arrow start point (slightly behind the position)
    arrow_start_x = arrow_x - 0.3 * arrow_length * dx_norm
    arrow_start_y = arrow_y - 0.3 * arrow_length * dy_norm
    
    # Arrow end point (slightly ahead of the position)
    arrow_end_x = arrow_x + 0.5 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.5 * arrow_length * dy_norm
    
    # Draw arrow using annotate
    ax.annotate('',  # Empty text (only arrow)
                xy=(arrow_end_x, arrow_end_y),   # Arrow head position
                xytext=(arrow_start_x, arrow_start_y),  # Arrow tail position
                arrowprops=dict(arrowstyle='->', 
                              color=f'C{i}',      # Different color for each arrow
                              linewidth=2.5,
                              mutation_scale=20,  # Controls arrow head size
                              shrinkA=0,          # No shrink at start
                              shrinkB=0),         # No shrink at end
                zorder=5)  # Ensure arrows are on top
    
    # Add a small dot at arrow position for clarity
    ax.scatter(arrow_x, arrow_y, color=f'C{i}', s=50, 
               edgecolor='black', linewidth=1, zorder=6, 
               label=f'Arrow {i+1} position')

# Add markers at the endpoints of the line
ax.scatter([x_start, x_end], [y_start, y_end], 
           c='red', s=100, zorder=5, 
           edgecolor='black', linewidth=2, 
           label='Endpoints')

# Add coordinate labels for endpoints
ax.text(x_start, y_start - 0.3, f'({x_start}, {y_start})', 
        ha='center', va='top', fontsize=10)
ax.text(x_end, y_end + 0.3, f'({x_end}, {y_end})', 
        ha='center', va='bottom', fontsize=10)

# Set axis properties
ax.set_xlim(-1, 6)
ax.set_ylim(-1, 4)
ax.set_aspect('equal')  # Keep aspect ratio 1:1
ax.grid(True, alpha=0.3, linestyle='--')

# Add title and labels
ax.set_title('Multiple Arrows on a Line Segment', fontsize=14, fontweight='bold')
ax.set_xlabel('X-axis', fontsize=12)
ax.set_ylabel('Y-axis', fontsize=12)

# Add legend
ax.legend(loc='upper left', fontsize=10)

# Add text explanation
info_text = f"Line length: {length:.2f} units\n"
info_text += f"Direction vector: ({dx}, {dy})\n"
info_text += f"Normalized direction: ({dx_norm:.2f}, {dy_norm:.2f})"
ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8))

plt.tight_layout()
plt.show()