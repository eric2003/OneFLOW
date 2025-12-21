import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

def draw_arrow_on_line(ax, x_start, y_start, x_end, y_end, position=0.5, 
                       arrow_style='->', color='blue', linewidth=2, 
                       head_size=15, zorder=2, show_connection=False):
    """
    Draw an arrow on a line segment at a specified position.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to draw on
    x_start, y_start : float
        Starting point coordinates
    x_end, y_end : float
        Ending point coordinates
    position : float (default=0.5)
        Position along the line where to place the arrow (0=start, 1=end)
    arrow_style : str (default='->')
        Arrow style (e.g., '->', '-|>', '<-', '<->', 'fancy')
    color : str or tuple (default='blue')
        Arrow color
    linewidth : float (default=2)
        Arrow line width
    head_size : float (default=15)
        Arrow head size (mutation_scale parameter)
    zorder : int (default=2)
        Drawing order (higher values are drawn on top)
    show_connection : bool (default=False)
        Whether to show the connection line (arrow shaft)
    """
    # Calculate the coordinates for the arrow position
    arrow_x = x_start + position * (x_end - x_start)
    arrow_y = y_start + position * (y_end - y_start)
    
    # Calculate direction vector
    dx = x_end - x_start
    dy = y_end - y_start
    
    # Normalize direction vector
    length = np.sqrt(dx**2 + dy**2)
    if length > 0:
        dx_norm = dx / length
        dy_norm = dy / length
    else:
        dx_norm, dy_norm = 0, 0
    
    # Define arrow length
    arrow_length = 0.2 * length
    
    # Calculate arrow start and end points
    arrow_start_x = arrow_x - 0.4 * arrow_length * dx_norm
    arrow_start_y = arrow_y - 0.4 * arrow_length * dy_norm
    arrow_end_x = arrow_x + 0.4 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.4 * arrow_length * dy_norm
    
    # Create arrow properties
    arrow_props = dict(arrowstyle=arrow_style,
                       color=color,
                       linewidth=linewidth,
                       mutation_scale=head_size,
                       shrinkA=0,
                       shrinkB=0)
    
    # If we don't want the connection line, set linestyle to 'none'
    if not show_connection:
        arrow_props['linestyle'] = 'none'
    
    # Draw the arrow
    arrow = ax.annotate('',
                        xy=(arrow_end_x, arrow_end_y),
                        xytext=(arrow_start_x, arrow_start_y),
                        arrowprops=arrow_props,
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
    
    draw_arrow_on_line(plt, 
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
