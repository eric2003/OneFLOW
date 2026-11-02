import matplotlib.pyplot as plt
import numpy as np

class BaseMesh:
    """Base class for mesh (main mesh/ghost mesh) with common grid logic"""
    def __init__(self, ncells, xstart, dx, lr=None):
        self.ncells = ncells
        self.nnodes = self.ncells + 1
        self.xstart = xstart  # Starting coordinate of the mesh
        self.dx = dx          # Grid spacing (negative value for left ghost mesh)
        self.lr = lr          # Identifier for ghost mesh: "L" (left) / "R" (right), None for main mesh
        
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
            self.y[i] = 0.0
        
        # Compute cell-center coordinates as midpoint of adjacent nodes
        for i in range(self.ncells):
            self.xcc[i] = 0.5 * (self.x[i] + self.x[i+1])
            self.ycc[i] = 0.0

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
        dy = 0.1 * abs(self.dx)  # Absolute value ensures positive vertical line length
        ipoints = [0,self.nnodes-1]
        for i in ipoints:
            xm = self.x[i]
            ym = self.y[i]
            plt.plot([xm, xm], [ym - 3*dy, ym + 3*dy], 'k-', linewidth=1)        

    def plot_vertical_interface_lines(self, indices=None):
        """Plot vertical lines at all mesh node positions (cell boundaries)"""
        dy = 0.1 * abs(self.dx)  # Absolute value ensures positive vertical line length
        if indices is None:
            indices = [i for i in range(0, self.nnodes)]
        for i in indices:
            xm = self.x[i]
            ym = self.y[i]
            plt.plot([xm, xm], [ym - dy, ym + dy], 'k-', linewidth=1)

class Ghost(BaseMesh):
    """Ghost mesh class with cell labeling and visualization"""
    def __init__(self, xstart, dx, ncells, lr):
        # Inherit initialization logic from BaseMesh class
        super().__init__(ncells=ncells, xstart=xstart, dx=dx, lr=lr)

    def plot(self):
        """Visualize ghost mesh: cell-center points and cell index labels"""
        ytext_shift = 0.5*abs(self.dx)  # Y-position for labels (avoid overlap with main mesh)
        for i in range(self.ncells):
            # Define cell labels based on left/right ghost mesh type
            if self.lr == "L":
                mystr = f"${- (i+1)}$"  # Label format for left ghost cells
            else:
                mystr = f"$N+{i+1}$"    # Label format for right ghost cells
            # Add text label at cell center
            plt.text(self.xcc[i], self.ycc[i]-ytext_shift, mystr, fontsize=12, ha='center')
        
        # Plot cell-center points (red fill with black edge)
        plt.scatter(self.xcc, self.ycc, s=50, facecolor='red', edgecolor='black', linewidth=1)
        
        indices = [i for i in range(1, self.nnodes)]
        self.plot_vertical_interface_lines(indices)
        
    def plot_ghost_vertical_interface_lines(self):
        """Plot vertical lines at all mesh node positions (cell boundaries)"""
        dy = 0.1 * abs(self.dx)  # Absolute value ensures positive vertical line length
        for i in range(1, self.nnodes-1):
            xm = self.x[i]
            ym = self.y[i]
            plt.plot([xm, xm], [ym - dy, ym + dy], 'k-', linewidth=1)        

class Mesh(BaseMesh):
    """Main mesh class with ghost mesh generation and management"""
    def __init__(self):
        # Define main mesh physical boundaries and grid resolution
        self.xmin = 0.0      # Left physical boundary of main mesh
        self.xmax = 1.0      # Right physical boundary of main mesh
        self.ncells = 9      # Number of cells in main mesh
        self.dx = (self.xmax - self.xmin) / self.ncells  # Grid spacing of main mesh
        self.nnodes = self.ncells + 1  # Number of nodes in main mesh
        
        # Initialize main mesh using BaseMesh constructor
        super().__init__(ncells=self.ncells, xstart=self.xmin, dx=self.dx, lr=None)
        self.nghosts = 2     # Number of ghost cell layers on each side

        # Create left ghost mesh (mirror extension to the left of main mesh)
        self.ghost_mesh_left = Ghost(
            xstart=self.xmin,
            dx=-self.dx,
            ncells=self.nghosts,
            lr="L"
        )
        
        # Create right ghost mesh (mirror extension to the right of main mesh)
        self.ghost_mesh_right = Ghost(
            xstart=self.xmax,
            dx=self.dx,
            ncells=self.nghosts,
            lr="R"
        )
        
    def generate_total_mesh(self):
        self.generate_mesh()
        self.ghost_mesh_left.generate_mesh()
        self.ghost_mesh_right.generate_mesh()

    def printinfo(self):
        """Print main mesh and ghost mesh information"""
        super().printinfo(prefix="Main Mesh")
        self.ghost_mesh_left.printinfo(prefix="Left Ghost Mesh")
        self.ghost_mesh_right.printinfo(prefix="Right Ghost Mesh")

    def plot(self):
        """Complete visualization of main mesh and ghost meshes"""
        # Plot main mesh cell-center points (black fill with black edge)
        plt.scatter(self.xcc, self.ycc, s=50, facecolor='black', edgecolor='black', linewidth=1)
        # Plot horizontal line connecting main mesh nodes
        plt.plot(self.x, self.y, 'k-', linewidth=1)
        # Plot vertical interface lines for main mesh
        indices = [i for i in range(1, self.nnodes-1)]
        self.plot_vertical_interface_lines(indices)
        self.plot_boundary_vertical_interface_lines()

        # Add boundary labels for main mesh
        dy = 0.1 * self.dx
        plt.text(self.x[0], self.y[0]+3*dy, r'$x=a$', fontsize=12, ha='center')
        plt.text(self.x[-1], self.y[0]+3*dy, r'$x=b$', fontsize=12, ha='center')
        
        # Plot ghost mesh components (points and labels)
        self.ghost_mesh_left.plot()
        self.ghost_mesh_right.plot()

def plot_cfd_figure():
    """Generate and save CFD mesh visualization figure"""
    # Configure LaTeX for text rendering and font settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif=['Times New Roman'])
    
    # Create figure with fixed dimensions
    plt.figure(figsize=(12, 4))

    # Initialize main mesh and generate grid coordinates
    mesh = Mesh()
    mesh.generate_total_mesh()
    
    # Print mesh information for verification
    mesh.printinfo()
    # Render all mesh components
    mesh.plot()

    # Set axis limits for symmetric display
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1, 1)
    
    # Set equal axis scale and hide axis lines
    plt.axis('equal')
    plt.axis('off')
    
    # Save figure with high resolution and tight bounding box
    plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
    # Display the figure
    plt.show()

if __name__ == '__main__':
    plot_cfd_figure()