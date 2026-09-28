from fractions import Fraction
import matplotlib.pyplot as plt
import numpy as np

def plot_vertical_boundary(border_x, y_start, y_end, lr):
    """
    Plot a vertical boundary with diagonal lines on one side.
   
    Parameters:
    -----------
    border_x : float
        The x-coordinate of the vertical boundary line
    y_start : float
        The starting y-coordinate of the vertical boundary
    y_end : float
        The ending y-coordinate of the vertical boundary
    lr : str
        Direction indicator: 'L' for left side, 'R' for right side
        Determines the orientation of diagonal lines
    """
    # Batch plot diagonal lines
    num_lines = 10 # Number of diagonal lines
   
    # Calculate the total height
    dy = y_end - y_start
   
    # Generate evenly spaced y positions for diagonal lines
    # Add small margins to avoid lines at the very edges
    y_positions = np.linspace(y_start + 0.02 * dy, y_end - 0.02 * dy, num_lines)
   
    # Length of each diagonal line (as percentage of total height)
    line_length = 0.2 * dy
   
    # Determine angle based on direction
    if lr == "L":
        angle = 180 + 30 # 210 degrees for left side
    else:
        angle = 30 # 30 degrees for right side
   
    # Convert angle to radians for trigonometric calculations
    angle_rad = np.deg2rad(angle)
   
    # Calculate dx and dy components for the diagonal lines
    dx = line_length * np.cos(angle_rad)
    dy_component = line_length * np.sin(angle_rad)
   
    # Plot the main vertical boundary line
    plt.plot([border_x, border_x], [y_start, y_end], 'k-', linewidth=2)
   
    # Plot each diagonal line
    for y in y_positions:
        # Start point: on the vertical line
        x1 = border_x
        y1 = y
       
        # End point: offset by dx and dy
        x2 = x1 + dx
        y2 = y1 + dy_component
       
        # Plot the diagonal line
        plt.plot([x1, x2], [y1, y2], color='blue', linewidth=1)

class BaseMesh:
    """Base class for mesh (main mesh/ghost mesh) with common grid logic"""
    def __init__(self, ncells, xstart, dx, ishift=0, ym=0, lr=None):
        self.ncells = ncells
        self.nnodes = self.ncells + 1
        self.xstart = xstart # Starting coordinate of the mesh
        self.dx = dx # Grid spacing (negative value for left ghost mesh)
        self.lr = lr # Identifier for ghost mesh: "L" (left) / "R" (right), None for main mesh
        self.ym = ym
        self.ishift = ishift
       
        # Initialize empty arrays for node coordinates and cell-center coordinates
        self.x = np.zeros(self.nnodes, dtype=np.float64)
        self.y = np.zeros(self.nnodes, dtype=np.float64)
        self.xcc = np.zeros(self.ncells, dtype=np.float64)
        self.ycc = np.zeros(self.ncells, dtype=np.float64)
    def generate_mesh(self):
        """Calculate node coordinates and cell-center coordinates for the mesh"""
        # Compute node coordinates along x-axis (y=0 for 1D mesh)
        for i in range(self.nnodes):
            self.x[i] = self.xstart + i * self.dx
            self.y[i] = self.ym
       
        # Compute cell-center coordinates as midpoint of adjacent nodes
        for i in range(self.ncells):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i+1])
            self.ycc[i] = 0.5 * (self.y[i] + self.y[i+1])
    def printinfo(self, prefix="Mesh"):
        """Print detailed mesh parameters and coordinates"""
        print(f"{prefix} ncells = {self.ncells}")
        print(f"{prefix} nnodes = {self.nnodes}")
        print(f"{prefix} xstart = {self.xstart:.6f}")
        if self.lr is not None:
            print(f"{prefix} lr = {self.lr}")
        print(f"{prefix} dx = {self.dx:.6f}")
        print(f"{prefix} x coordinates = {self.x}")
        print(f"{prefix} cell-center x coordinates = {self.xcc}")
       
    def plot_boundary_vertical_interface_lines(self):
        """Plot vertical lines at all mesh node positions (cell boundaries)"""
        dy = 0.1 * abs(self.dx) # Absolute value ensures positive vertical line length
        ipoints = [0,self.nnodes-1]
        for i in ipoints:
            xm = self.x[i]
            ym = self.y[i]
            #plt.plot([xm, xm], [ym - 3*dy, ym + 3*dy], 'k-', linewidth=1)
            lr = "L" if i == 0 else "R"
            plot_vertical_boundary(xm,ym - 3*dy,ym + 3*dy, lr)
    def plot_vertical_interface_lines(self, indices=None):
        """Plot vertical lines at all mesh node positions (cell boundaries)"""
        dy = 0.1 * abs(self.dx) # Absolute value ensures positive vertical line length
        if indices is None:
            indices = [i for i in range(0, self.nnodes)]
        for i in indices:
            xm = self.x[i]
            ym = self.y[i]
            plt.plot([xm, xm], [ym - dy, ym + dy], 'k-', linewidth=1)

class Ghost(BaseMesh):
    """Ghost mesh class with cell labeling and visualization"""
    def __init__(self, xstart, dx, ncells, ishift, ym, lr):
        # Inherit initialization logic from BaseMesh class
        super().__init__(ncells=ncells, xstart=xstart, dx=dx, ishift=ishift, ym=ym, lr=lr)
       
    def plot_cell_label(self):
        ytext_shift = 0.5*abs(self.dx) # Y-position for labels (avoid overlap with main mesh)
        print(f"self.ishift={self.ishift}")
        for i in range(self.ncells):
            # Define cell labels based on left/right ghost mesh type
            if self.lr == "L":
                cell_label = f"${- i+self.ishift}$" # Label format for left ghost cells
            else:
                ii = i+1+self.ishift
                if ii == 0:
                    cell_label = f"$N$"
                else:
                    cell_label = f"$N+{i+1+self.ishift}$" # Label format for right ghost cells
            # Add text label at cell center
            plt.text(self.xcc[i], self.ycc[i]-ytext_shift, cell_label, fontsize=12, ha='center')
           
    def plot(self):
        """Visualize ghost mesh: cell-center points and cell index labels"""
        self.plot_cell_label()
        # Plot cell-center points (red fill with black edge)
        plt.scatter(self.xcc, self.ycc, s=50, facecolor='red', edgecolor='black', linewidth=1)
       
        indices = [i for i in range(1, self.nnodes)]
        self.plot_vertical_interface_lines(indices)

class Mesh(BaseMesh):
    """Main mesh class with ghost mesh generation and management"""
    def __init__(self, nghosts=2,ishift=0, ym=0, periodic=False): # Added: periodic parameter, default False
        print(f"Mesh ishift={ishift}")
        self.ym = ym
        self.periodic = periodic # Added: periodic flag
        self.xmin = 0.0 # Left physical boundary of main mesh
        self.xmax = 1.0 # Right physical boundary of main mesh
        #self.ncells = 9 # Number of cells in main mesh
        self.ncells = 7 # Number of cells in main mesh
        self.dx = (self.xmax - self.xmin) / self.ncells # Grid spacing of main mesh
        self.nnodes = self.ncells + 1 # Number of nodes in main mesh
       
        # Initialize main mesh using BaseMesh constructor
        super().__init__(ncells=self.ncells, xstart=self.xmin, dx=self.dx, ishift=ishift, ym=self.ym, lr=None)
        print(f"Mesh self.ishift={self.ishift}")
        #self.nghosts = 2 # Number of ghost cell layers on each side
        self.nghosts = nghosts # Number of ghost cell layers on each side
        # Create left ghost mesh (mirror extension to the left of main mesh)
        self.ghost_mesh_left = Ghost(
            xstart=self.xmin,
            dx=-self.dx,
            ncells=self.nghosts,
            ishift=ishift,
            ym=self.ym,
            lr="L"
        )
       
        # Create right ghost mesh (mirror extension to the right of main mesh)
        self.ghost_mesh_right = Ghost(
            xstart=self.xmax,
            dx=self.dx,
            ncells=self.nghosts,
            ishift=ishift,
            ym=self.ym,
            lr="R"
        )
       
        self.generate_total_mesh()
   
        # Print mesh information for verification
        self.printinfo()
        # Render all mesh components
        self.plot()
       
    def generate_total_mesh(self):
        self.generate_mesh()
        self.ghost_mesh_left.generate_mesh()
        self.ghost_mesh_right.generate_mesh()
    def printinfo(self):
        """Print main mesh and ghost mesh information"""
        super().printinfo(prefix="Main Mesh")
        self.ghost_mesh_left.printinfo(prefix="Left Ghost Mesh")
        self.ghost_mesh_right.printinfo(prefix="Right Ghost Mesh")
       
    def getstringFrac(self, a):
        absa = abs(a)
        sign = "-" if a < 0 else ""
        str_a = f"{{{absa.numerator}}}{{{absa.denominator}}}"
        str_a = rf"$x_{{{sign}\frac{str_a}}}$"
        return str_a
       
    def getstringNFrac(self, a):
        absa = abs(a)
        sign = "-" if a < 0 else "+"
        str_na = f"{{{absa.numerator}}}{{{absa.denominator}}}"
        str_na = rf"$x_{{N{sign}\frac{str_na}}}$"
        return str_na
       
    def get_label(self,strin,index):
        absindex = abs(index)
        sign = "-" if index < 0 else "+"
        if index == 0:
            #label = f"${strin}$"
            label = f"{strin}"
        else:
            #label = f"${strin}{sign}{absindex}$"
            label = f"{strin}{sign}{absindex}"
        return label
       
    def get_dollar_label(self,strin,index):
        return f"${self.get_label(strin,index)}$"
    def plot_cell_label(self):
        dytext = 0.5*abs(self.dx)
        self.nlabels = 2
        for i in range(self.ncells):
            if i < self.nlabels:
                  cell_label = f"${i+1+self.ishift}$"
            elif i > self.ncells - 1 - self.nlabels:
                inew = i - (self.ncells - 1) + self.ishift
                cell_label = self.get_dollar_label("N",inew)
            else:
                cell_label=""
            # Add text label at cell center
            plt.text(self.xcc[i], self.ycc[i]-dytext, cell_label, fontsize=12, ha='center')
        icenter = self.ncells // 2
        plt.text(self.xcc[icenter], self.ycc[icenter]-dytext, f"$i$", fontsize=12, ha='center')
        cell_label = self.get_label("N",self.ishift+1)
        #text = f"$ied-ist=N$"
        #plt.text(self.xcc[icenter], self.ycc[icenter]-3*dytext, text, fontsize=12, ha='center')
        xm = self.xcc[self.ncells-1] + self.dx
        ym = self.ycc[self.ncells-1]
        #plt.arrow(self.xcc[0], self.ycc[0]-2.5*dytext, 0, 1.0*dytext, head_width=0.01, head_length=0.01, fc='blue', ec='blue', linewidth=1)
        #plt.arrow(xm, ym-2.5*dytext, 0, 1.0*dytext, head_width=0.01, head_length=0.01, fc='blue', ec='blue', linewidth=1)
        text1 = f"$ist={1+self.ishift}$"
        text2 = f"$ied={cell_label}$"
        #plt.text(self.xcc[0], self.ycc[0]-3*dytext, text1, fontsize=12, ha='center')
        #plt.text(xm, ym-3*dytext, text2, fontsize=12, ha='center')
    def plot_cell_center(self):
        # Plot main mesh cell-center points (black fill with black edge)
        xcc_new = []
        ycc_new = []
        for i in range(self.ncells):
            if self.cellmark[i] == 1:
                xcc_new.append( self.xcc[i] )
                ycc_new.append( self.ycc[i] )
        plt.scatter(xcc_new, ycc_new, s=50, facecolor='black', edgecolor='black', linewidth=1)
       
    def plot_cell_mesh(self):
        self.icenter = self.ncells // 2
        self.nlabels = 2
        self.cellmark = np.zeros(self.ncells, dtype=int)
        for i in range(self.ncells):
            if i < self.nlabels:
                self.cellmark[i] = 1
            elif i > self.ncells - 1 - self.nlabels:
                self.cellmark[i] = 1
        self.cellmark[self.icenter] = 1
       
        self.plot_cell_center()
        self.plot_cell_label()
       
        # Plot horizontal line connecting main mesh nodes
        for i in range(self.ncells):
            if self.cellmark[i] == 1:
                plt.plot([self.x[i], self.x[i+1]], [self.y[i], self.y[i+1]], 'k-', linewidth=1)
            else:
                plt.plot([self.x[i], self.x[i+1]], [self.y[i], self.y[i+1]], 'k--', linewidth=1)
               
       
    def plot(self):
        """Complete visualization of main mesh and ghost meshes"""
       
        self.plot_cell_mesh()
        # Plot vertical interface lines for main mesh
        indices = [i for i in range(1, self.nnodes-1)]
        self.plot_vertical_interface_lines(indices)
        self.plot_boundary_vertical_interface_lines()
       
        # Add boundary labels for main mesh
        dym = abs(self.dx)
        plt.text(self.x[0], self.y[0]+0.5*dym, r'$x=a$', fontsize=12, ha='center')
        plt.text(self.x[-1], self.y[0]+0.5*dym, r'$x=b$', fontsize=12, ha='center')
       
        half = Fraction(1,2)
        a1 = half+self.ishift
        a2 = half+self.ishift+1
        print(f"a1={a1}")
        print(f"a2={a2}")
        str_a1 = self.getstringFrac(a1)
        str_a2 = self.getstringFrac(a2)
       
        an1 = half+self.ishift
        an2 = -half+self.ishift
       
        str_na1 = self.getstringNFrac(an1)
        str_na2 = self.getstringNFrac(an2)
       
        print(f"an1={an1}")
        print(f"an2={an2}")
       
        print(f"str_na1={str_na1}")
        print(f"str_na2={str_na2}")
       
        #plt.text(self.x[0], self.y[0]-dym, str_a1, fontsize=12, ha='center')
        #plt.text(self.x[1], self.y[1]-dym, str_a2, fontsize=12, ha='center')
        #plt.text(self.x[-1], self.y[-1]-dym, str_na1, fontsize=12, ha='center')
        #plt.text(self.x[-2], self.y[-2]-dym, str_na2, fontsize=12, ha='center')
       
        # Plot ghost mesh components (points and labels)
        self.ghost_mesh_left.plot()
        self.ghost_mesh_right.plot()
        # Added: If periodic=True, draw periodic connections
        if self.periodic:
            self.plot_periodic_connections()
       
    # Added method: Plot periodic boundary fold-line arrows from source cells to ghost cells, indicating value assignment
    def plot_periodic_connections(self):
        """Plot periodic boundary fold-line arrows from source cells to ghost cells, indicating value assignment"""
        #vlen = 0.3 * abs(self.dx)  # Vertical line length
        vlen = 0.6 * abs(self.dx)  # Vertical line length
        offset_y = 0.05 * abs(self.dx)  # Small offset from cell center to avoid overlap
        #h_offset = 0.02  # Horizontal line fine-tuning to avoid overlap
        h_offset = 0.0
        head_width = 0.01  # Arrow head width
        head_length = 0.01  # Arrow head length
        lw = 1  # Line width

        # Left ghost assignments: Draw above x-axis (y > 0)
        # cell N (main i=6) -> u[0] (left ghost i=1), cell N-1 (main i=5) -> u[1] (left ghost i=0)
        connections_left = [(6, 0, '-', 'blue'), (5, 1, '--', 'purple')]  # Different line styles and colors
        for idx, (src_i, tgt_i, ls, color) in enumerate(connections_left):
            #print(f"连接 {idx}: 从节点{src_i}到节点{tgt_i}, 线型{ls}, 颜色{color}")
    
            src_x, src_y = self.xcc[src_i], self.ycc[src_i]
            tgt_x, tgt_y = self.ghost_mesh_left.xcc[tgt_i], self.ghost_mesh_left.ycc[tgt_i]
           
            # Source vertical line upward with arrow (from src_y + offset_y to src_y + vlen + offset_y)
            dy = idx * vlen
            start_y_src = src_y + offset_y
            end_y_src = src_y + vlen + offset_y + dy
            plt.arrow(src_x, start_y_src, 0, end_y_src - start_y_src, 
                      head_width=head_width, head_length=head_length, fc=color, ec=color, linewidth=lw)
            # Horizontal line to target above
            #plt.plot([src_x, tgt_x], [end_y_src, tgt_y + vlen + offset_y], color=color, ls=ls, linewidth=lw)
            plt.plot([src_x, tgt_x], [end_y_src, end_y_src], color=color, ls=ls, linewidth=lw)
            print(f"end_y_src,tgt_y + vlen + offset_y={end_y_src,end_y_src}")
            # Target vertical line downward with arrow (from tgt_y + vlen + offset_y to tgt_y + offset_y)
            start_y_tgt = tgt_y + vlen + offset_y + dy
            #start_y_tgt = start_y_src
            end_y_tgt = tgt_y + offset_y
            #print(f"")
            plt.arrow(tgt_x, start_y_tgt, 0, end_y_tgt - start_y_tgt, 
                      head_width=head_width, head_length=head_length, fc=color, ec=color, linewidth=lw)

        # Right ghost assignments: Draw below x-axis (y < 0) by mirroring the logic
        # cell2 (main i=1) -> N+2 (right ghost i=0), cell3 (main i=2) -> N+3 (right ghost i=1)
        connections_right = [(0, 0, '-', 'green'), (1, 1, '--', 'orange')]  # Different line styles and colors
        for src_i, tgt_i, ls, color in connections_right:
            src_x, src_y = self.xcc[src_i], self.ycc[src_i]
            tgt_x, tgt_y = self.ghost_mesh_right.xcc[tgt_i], self.ghost_mesh_right.ycc[tgt_i]
           
            # Source vertical line downward with arrow (from src_y - offset_y to src_y - vlen - offset_y)
            start_y_src = src_y - offset_y
            end_y_src = src_y - vlen - offset_y
            plt.arrow(src_x, start_y_src, 0, end_y_src - start_y_src, 
                      head_width=head_width, head_length=head_length, fc=color, ec=color, linewidth=lw)
            # Horizontal line to target below
            plt.plot([src_x, tgt_x], [end_y_src, end_y_src], color=color, ls=ls, linewidth=lw)
            # Target vertical line upward with arrow (from tgt_y - vlen - offset_y to tgt_y - offset_y)
            start_y_tgt = tgt_y - vlen - offset_y + h_offset
            end_y_tgt = tgt_y - offset_y
            plt.arrow(tgt_x, start_y_tgt, 0, end_y_tgt - start_y_tgt, 
                      head_width=head_width, head_length=head_length, fc=color, ec=color, linewidth=lw)

def plot_cfd_figure(periodic=False): # Added: periodic parameter, default False
    """Generate and save CFD mesh visualization figure"""
    # Configure LaTeX for text rendering and font settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif=['Times New Roman'])
   
    # Create figure with fixed dimensions
    plt.figure(figsize=(12, 8))
    # Initialize main mesh and generate grid coordinates
    nghost2 = 2
    mesh = Mesh(nghost2,nghost2-1, periodic=periodic) # Pass periodic
   
    # Set axis limits for symmetric display
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1, 1)
   
    # Set equal axis scale and hide axis lines
    plt.axis('equal')
    plt.axis('off')
   
    # Save figure with high resolution and tight bounding box
    filename = 'cfd_periodic.png' if periodic else 'cfd.png'
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    # Display the figure
    plt.show()

if __name__ == '__main__':
    plot_cfd_figure(periodic=True) # Example: Enable periodic visualization