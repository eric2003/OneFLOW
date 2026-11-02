import matplotlib.pyplot as plt
import numpy as np

def draw_arrow_on_line(ax, x_start, y_start, x_end, y_end, position=0.5, 
                       arrow_style='->', color='blue', linewidth=2, 
                       head_size=15, label=None, zorder=2):
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
    label : str or None (default=None)
        Label for the arrow (for legend)
    zorder : int (default=2)
        Drawing order (higher values are drawn on top)
    
    Returns:
    --------
    arrow_annotation : matplotlib.text.Annotation
        The arrow annotation object
    """
    # Validate position parameter
    if position < 0 or position > 1:
        raise ValueError("Position must be between 0 and 1")
    
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
    
    # Define arrow length (relative to line length)
    arrow_length = 0.2 * length  # 20% of line length
    
    # Calculate arrow start and end points
    # Arrow points in the direction of the line
    arrow_start_x = arrow_x - 0.4 * arrow_length * dx_norm
    arrow_start_y = arrow_y - 0.4 * arrow_length * dy_norm
    arrow_end_x = arrow_x + 0.4 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.4 * arrow_length * dy_norm
    
    # Draw the arrow
    arrow = ax.annotate('',
                        xy=(arrow_end_x, arrow_end_y),
                        xytext=(arrow_start_x, arrow_start_y),
                        arrowprops=dict(arrowstyle=arrow_style,
                                       color=color,
                                       linewidth=linewidth,
                                       mutation_scale=head_size,
                                       shrinkA=0,
                                       shrinkB=0),
                        zorder=zorder)
    
    # Add a small marker at the arrow position
    marker = ax.scatter(arrow_x, arrow_y, color=color, s=30, 
                       zorder=zorder+1, alpha=0.6)
    
    # Add label if provided
    if label:
        ax.text(arrow_x, arrow_y + 0.05*length, label, 
                color=color, fontsize=9, ha='center', va='bottom',
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
    
    return arrow, marker

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
