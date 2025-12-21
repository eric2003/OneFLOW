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
    arrow_end_x = arrow_x + 0.6 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.6 * arrow_length * dy_norm
    
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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Example 1: Single line with arrows at different positions
    ax1.set_title("Arrows at Different Positions", fontsize=14, fontweight='bold')
    
    # Draw the line
    line_x = [0, 5]
    line_y = [0, 3]
    ax1.plot(line_x, line_y, 'gray', linewidth=3, alpha=0.3, label='Base line')
    
    # Draw arrows at different positions
    positions = [0.0, 0.25, 0.5, 0.75, 1.0]
    colors = ['red', 'orange', 'green', 'blue', 'purple']
    labels = ['Start (0.0)', '0.25', 'Middle (0.5)', '0.75', 'End (1.0)']
    
    for pos, color, label in zip(positions, colors, labels):
        draw_arrow_on_line(ax1, 
                          x_start=line_x[0], y_start=line_y[0],
                          x_end=line_x[1], y_end=line_y[1],
                          position=pos,
                          arrow_style='->',
                          color=color,
                          linewidth=2,
                          head_size=15,
                          label=label,
                          zorder=5)
    
    # Mark endpoints
    ax1.scatter(line_x, line_y, c='black', s=100, zorder=10, label='Endpoints')
    
    # Add coordinate labels
    ax1.text(line_x[0], line_y[0]-0.4, f'({line_x[0]}, {line_y[0]})', 
            ha='center', fontsize=10)
    ax1.text(line_x[1], line_y[1]+0.4, f'({line_x[1]}, {line_y[1]})', 
            ha='center', fontsize=10)
    
    ax1.set_xlim(-1, 6)
    ax1.set_ylim(-1, 5)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left')
    ax1.set_xlabel('X-axis')
    ax1.set_ylabel('Y-axis')
    
    # Example 2: Multiple lines with arrows
    ax2.set_title("Multiple Lines with Arrows", fontsize=14, fontweight='bold')
    
    # Define multiple lines
    lines = [
        {'start': (0, 0), 'end': (4, 1), 'color': 'blue', 'label': 'Line 1'},
        {'start': (1, 0), 'end': (3, 3), 'color': 'green', 'label': 'Line 2'},
        {'start': (0, 2), 'end': (5, 0), 'color': 'red', 'label': 'Line 3'},
    ]
    
    # Draw each line with arrows
    for i, line in enumerate(lines):
        x1, y1 = line['start']
        x2, y2 = line['end']
        color = line['color']
        label = line['label']
        
        # Draw the line
        ax2.plot([x1, x2], [y1, y2], color=color, linewidth=3, 
                alpha=0.3, label=label)
        
        # Draw arrows at three positions
        for pos in [0.2, 0.5, 0.8]:
            draw_arrow_on_line(ax2,
                              x_start=x1, y_start=y1,
                              x_end=x2, y_end=y2,
                              position=pos,
                              arrow_style='-|>',
                              color=color,
                              linewidth=1.5,
                              head_size=12,
                              zorder=5)
        
        # Mark endpoints
        ax2.scatter([x1, x2], [y1, y2], c=color, s=80, zorder=10, edgecolors='black')
    
    ax2.set_xlim(-0.5, 5.5)
    ax2.set_ylim(-0.5, 3.5)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    ax2.set_xlabel('X-axis')
    ax2.set_ylabel('Y-axis')
    
    plt.tight_layout()
    plt.show()

# Example 3: Interactive demonstration
def interactive_demo():
    """Interactive demonstration with user-specified position."""
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Define a line
    x_start, y_start = 1, 1
    x_end, y_end = 8, 6
    
    # Draw the base line
    ax.plot([x_start, x_end], [y_start, y_end], 
            'gray', linewidth=4, alpha=0.2, label='Base line')
    
    # Test different arrow styles at position 0.5
    arrow_styles = ['->', '-|>', '<-', '<->', 'fancy', 'simple']
    colors = plt.cm.Set2(np.linspace(0, 1, len(arrow_styles)))
    
    for i, (style, color) in enumerate(zip(arrow_styles, colors)):
        # Offset each arrow slightly for visibility
        offset = 0.1 * i
        position = 0.3 + offset
        
        if position <= 1.0:  # Ensure position is valid
            arrow, marker = draw_arrow_on_line(ax,
                                             x_start=x_start, y_start=y_start,
                                             x_end=x_end, y_end=y_end,
                                             position=position,
                                             arrow_style=style,
                                             color=color,
                                             linewidth=2,
                                             head_size=15,
                                             label=f"style='{style}'",
                                             zorder=5)
    
    # Mark endpoints
    ax.scatter([x_start, x_end], [y_start, y_end], 
               c='black', s=150, zorder=10, label='Endpoints')
    
    # Add coordinate labels
    ax.text(x_start, y_start-0.5, f'Start: ({x_start}, {y_start})', 
            ha='center', fontsize=11, fontweight='bold')
    ax.text(x_end, y_end+0.5, f'End: ({x_end}, {y_end})', 
            ha='center', fontsize=11, fontweight='bold')
    
    # Calculate and display line information
    dx = x_end - x_start
    dy = y_end - y_start
    length = np.sqrt(dx**2 + dy**2)
    
    info_text = f"Line Information:\n"
    info_text += f"• Length: {length:.2f} units\n"
    info_text += f"• Direction: ({dx:.1f}, {dy:.1f})\n"
    info_text += f"• Angle: {np.degrees(np.arctan2(dy, dx)):.1f}°"
    
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9))
    
    # Set plot properties
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_title("Arrow Styles at Different Positions", fontsize=16, fontweight='bold')
    ax.set_xlabel('X-axis', fontsize=12)
    ax.set_ylabel('Y-axis', fontsize=12)
    ax.legend(loc='lower right', fontsize=10)
    
    plt.tight_layout()
    plt.show()

# Run examples
if __name__ == "__main__":
    print("Example 1: Arrows at different positions")
    example_usage()
    
    print("\nExample 2: Different arrow styles")
    interactive_demo()