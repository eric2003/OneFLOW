import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

def draw_arrow_only(ax, x_start, y_start, x_end, y_end, position=0.5,
                    arrow_style='->', color='blue', linewidth=2,
                    head_size=15, zorder=2):
    """
    Draw only the arrow head without the connecting line.
    """
    # Calculate arrow position
    arrow_x = x_start + position * (x_end - x_start)
    arrow_y = y_start + position * (y_end - y_start)
    
    # Calculate direction
    dx = x_end - x_start
    dy = y_end - y_start
    length = np.sqrt(dx**2 + dy**2)
    
    if length > 0:
        dx_norm = dx / length
        dy_norm = dy / length
    else:
        dx_norm, dy_norm = 0, 0
    
    # Very small offset to create just the arrow head
    offset = 0.001 * length
    
    # Create arrow
    arrow = ax.annotate('',
                        xy=(arrow_x + offset * dx_norm, arrow_y + offset * dy_norm),
                        xytext=(arrow_x - offset * dx_norm, arrow_y - offset * dy_norm),
                        arrowprops=dict(arrowstyle=arrow_style,
                                       color=color,
                                       linewidth=linewidth,
                                       mutation_scale=head_size,
                                       #linestyle='none',  # No line!
                                       shrinkA=0,
                                       shrinkB=0),
                        zorder=zorder)
    
    return arrow
# Example usage
def example_usage():
    """Example demonstrating the draw_arrow_on_line function."""
    
    # Create figure
    plt.figure(figsize=(8, 8))
    
    # Draw the line
    line_x = [0, 1]
    line_y = [0, 1]
    plt.plot(line_x, line_y, 'gray', linewidth=3, alpha=0.3, label='Base line')
    
    draw_arrow_only(plt, 
                       x_start=line_x[0], y_start=line_y[0],
                       x_end=line_x[1], y_end=line_y[1]
                      )
    
    
    plt.xlim(-0.5, 5.5)
    plt.ylim(-0.5, 3.5)
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    example_usage()
