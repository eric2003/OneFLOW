import numpy as np
import matplotlib.pyplot as plt

def read_results(filename):
    """Read results from Fortran output file"""
    try:
        data = np.loadtxt(filename, comments='#')
        return data[:, 0], data[:, 1]
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return np.array([]), np.array([])

def main():
    print("Reading Fortran output files...")
    
    # Read all data files
    x_eno, u_eno = read_results('eno_results.txt')
    x_weno, u_weno = read_results('weno_results.txt')
    x_analytical, u_analytical = read_results('analytical_results.txt')
    
    # Check if we have data
    if len(x_eno) == 0 or len(x_weno) == 0 or len(x_analytical) == 0:
        print("Error: Could not read all data files.")
        print("Make sure to run the Fortran program first.")
        return
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    # Plot results
    plt.plot(x_eno, u_eno, 'bo-', linewidth=1, markersize=3, 
             markerfacecolor='none', label='ENO3')
    plt.plot(x_weno, u_weno, 'gs-', linewidth=1, markersize=3,
             markerfacecolor='none', label='WENO3')
    plt.plot(x_analytical, u_analytical, 'r-', linewidth=2, label='Analytical')
    
    # Customize plot
    plt.title('1D Convection: ENO3 vs WENO3 (t=0.625)')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save and show
    plt.tight_layout()
    plt.savefig('comparison.png', dpi=150)
    plt.show()
    
    print("Plot saved as comparison.png")

if __name__ == "__main__":
    main()