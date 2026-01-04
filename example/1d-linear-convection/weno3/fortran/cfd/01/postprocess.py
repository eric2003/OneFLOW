import numpy as np
import matplotlib.pyplot as plt

def read_results(filename):
    """Read results from Fortran output file"""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1]

def main():
    # Read results
    x_eno, u_eno = read_results('eno_results.txt')
    x_weno, u_weno = read_results('weno_results.txt')
    x_analytical, u_analytical = read_results('analytical_results.txt')
    
    # Create plot
    plt.figure(figsize=(12, 8))
    plt.plot(x_eno, u_eno, 'bo-', linewidth=1.5, markersize=4, 
             markerfacecolor='none', label='ENO3 (Rusanov)')
    plt.plot(x_weno, u_weno, 'gs-', linewidth=1.5, markersize=4,
             markerfacecolor='none', label='WENO3 (Rusanov)')
    plt.plot(x_analytical, u_analytical, 'r--', linewidth=2, label='Analytical')
    
    plt.title('1D Convection Equation: ENO3 vs WENO3 Comparison (t=0.625)', fontsize=14)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('u', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save figure
    plt.savefig('eno_weno_comparison.png', dpi=300, bbox_inches='tight')
    print("Plot saved as: eno_weno_comparison.png")
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    main()