import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve

def solve_convection_diffusion(a=0.0, nu=0.0, t0=0.1, total_time=1.0, 
                              xmin=-5.0, xmax=10.0, Nx=200, Nt=1000):
    """
    Solve 1D convection-diffusion equation (supports toggling convection and diffusion terms)
    Equation: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    
    Parameters:
    a: convection velocity (set to 0 to turn off convection)
    nu: diffusion coefficient (set to 0 to turn off diffusion)
    t0: width parameter for initial Gaussian distribution
    total_time: total simulation time
    xmin, xmax: spatial domain boundaries
    Nx: number of spatial grid points
    Nt: number of time steps
    """
    dx = (xmax - xmin) / Nx
    dt = total_time / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    
    # Initial condition: Gaussian distribution
    if nu > 0:
        u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    else:
        # For pure convection (nu=0), use a small regularization
        u = np.exp(-x**2 / (4*0.01*t0)) / np.sqrt(4*np.pi*0.01*t0)
    
    # Time stepping (Crank-Nicolson scheme)
    for n in range(Nt):
        N = Nx + 1  # total number of points
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        b = np.zeros(N)
        
        # Interior points (1 to Nx-1)
        for i in range(1, Nx):
            # Convection coefficient (central difference)
            conv_coef = a / (2*dx) if a != 0 else 0.0
            # Diffusion coefficient (central second-order difference)
            diff_coef = nu / (dx**2) if nu != 0 else 0.0
            
            # Implicit matrix coefficients (I + 0.5*dt*L)
            main_diag[i] = 1.0 + dt * diff_coef
            lower_diag[i] = -0.5 * dt * (conv_coef + diff_coef)
            upper_diag[i] = 0.5 * dt * (conv_coef - diff_coef)
            
            # Right-hand side vector ((I - 0.5*dt*L)u^n)
            b[i] = (1.0 - dt * diff_coef) * u[i] + \
                   0.5 * dt * (conv_coef + diff_coef) * u[i-1] + \
                   0.5 * dt * (-conv_coef + diff_coef) * u[i+1]
        
        # Boundary conditions: zero gradient (Neumann)
        main_diag[0] = 1.0
        upper_diag[0] = -1.0
        b[0] = 0.0
        
        main_diag[Nx] = 1.0
        lower_diag[Nx] = -1.0
        b[Nx] = 0.0
        
        # Build and solve sparse matrix
        A = sparse.diags([lower_diag[1:], main_diag, upper_diag[:-1]], 
                        [-1, 0, 1], format='csr')
        u = spsolve(A, b)
    
    return x, u

def exact_solution(x, t, a=0.0, nu=0.0, t0=0.1):
    """
    Exact solution of convection-diffusion equation
    u(x,t) = 1/√[4πν(t+t0)] * exp[-(x-at)²/(4ν(t+t0))]
    For pure convection (nu=0), the Gaussian becomes a delta function
    For numerical stability, we use a small nu when nu=0
    """
    if nu == 0:
        # For pure convection, use a small diffusion for visualization
        nu_small = 0.01
        return 1.0 / np.sqrt(4*np.pi*nu_small*(t+t0)) * np.exp(-(x-a*t)**2/(4*nu_small*(t+t0)))
    else:
        return 1.0 / np.sqrt(4*np.pi*nu*(t+t0)) * np.exp(-(x-a*t)**2/(4*nu*(t+t0)))

def plot_comparison(cases, total_time=1.0):
    """
    Plot comparison for multiple parameter sets
    
    Parameters:
    cases: list, each element is (a, nu, label, color, linestyle)
    total_time: simulation time
    """
    plt.figure(figsize=(15, 8))
    
    for i, (a, nu, label, color, ls) in enumerate(cases):
        # Numerical solution
        x, u_num = solve_convection_diffusion(a=a, nu=nu, total_time=total_time)
        
        # Exact solution
        u_exact = exact_solution(x, t=total_time, a=a, nu=nu)
        
        # Plot comparison
        plt.subplot(2, 3, i+1)
        plt.plot(x, u_num, color=color, linestyle='-', linewidth=2, label='Numerical')
        plt.plot(x, u_exact, color=color, linestyle=ls, linewidth=2, alpha=0.7, label='Exact')
        
        plt.xlabel('x')
        plt.ylabel('u(x,T)')
        plt.title(f'{label}\na={a}, ν={nu}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Calculate and display error
        error = np.max(np.abs(u_num - u_exact))
        plt.text(0.02, 0.95, f'Max error: {error:.2e}', 
                transform=plt.gca().transAxes, fontsize=9,
                verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.show()

def main():
    """
    Main function: test convection-diffusion effects under different cases
    """
    print("="*60)
    print("Convection-Diffusion Equation Test: Gaussian Initial Distribution")
    print("Equation: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²")
    print("="*60)
    
    # Define test cases
    cases = [
        # (a, nu, label, color, linestyle)
        (0.0, 0.1, "Pure Diffusion", "blue", "--"),      # diffusion only
        (1.0, 0.01, "Pure Convection", "red", "--"),     # convection only (small nu for visualization)
        (1.0, 0.1, "Convection+Diffusion", "green", "--"),  # both terms
        (2.0, 0.05, "Strong Convection", "purple", "--"),
        (0.5, 0.2, "Strong Diffusion", "orange", "--"),
        (1.0, 0.01, "Convection Dominated", "brown", "--"),
    ]
    
    # Plot comparison
    plot_comparison(cases)
    
    # Detailed analysis for each case
    print("\nDetailed Analysis:")
    for a, nu, label, color, ls in cases:
        print(f"\n--- {label} (a={a}, ν={nu}) ---")
        
        # Compute numerical and exact solutions
        x, u_num = solve_convection_diffusion(a=a, nu=nu)
        u_exact = exact_solution(x, t=1.0, a=a, nu=nu)
        
        # Calculate errors
        error = np.abs(u_num - u_exact)
        max_error = np.max(error)
        rms_error = np.sqrt(np.mean(error**2))
        
        print(f"Maximum absolute error: {max_error:.2e}")
        print(f"RMS error: {rms_error:.2e}")
        
        # Plot individual comparison
        plt.figure(figsize=(10, 6))
        
        plt.subplot(2, 1, 1)
        plt.plot(x, u_num, color=color, linestyle='-', linewidth=2, label='Numerical')
        plt.plot(x, u_exact, color=color, linestyle=ls, linewidth=2, alpha=0.7, label='Exact')
        plt.xlabel('x')
        plt.ylabel('u(x,T)')
        plt.title(f'{label}: a={a}, ν={nu}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(2, 1, 2)
        plt.plot(x, error, 'r-', linewidth=1.5)
        plt.xlabel('x')
        plt.ylabel('Absolute Error')
        plt.title(f'Error Distribution (Max error: {max_error:.2e})')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()

def quick_test():
    """Quick test of three basic cases"""
    plt.figure(figsize=(15, 4))
    
    test_cases = [
        ("Pure Diffusion", 0.0, 0.1, "blue"),
        ("Pure Convection", 1.0, 0.01, "red"),
        ("Convection+Diffusion", 1.0, 0.1, "green"),
    ]
    
    for i, (title, a, nu, color) in enumerate(test_cases, 1):
        # Compute numerical solution
        x, u_num = solve_convection_diffusion(a=a, nu=nu)
        
        # Compute exact solution
        u_exact = exact_solution(x, t=1.0, a=a, nu=nu)
        
        # Plot
        plt.subplot(1, 3, i)
        plt.plot(x, u_num, color=color, linestyle='-', linewidth=2, label='Numerical')
        plt.plot(x, u_exact, color=color, linestyle='--', linewidth=2, alpha=0.7, label='Exact')
        
        plt.xlabel('x')
        plt.ylabel('u(x,T)')
        plt.title(title)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Calculate and display error
        error = np.max(np.abs(u_num - u_exact))
        plt.text(0.02, 0.95, f'Max error: {error:.2e}', 
                transform=plt.gca().transAxes, fontsize=10,
                verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Select run mode:")
    print("1. Full test (6 cases)")
    print("2. Quick test (3 basic cases)")
    
    try:
        choice = input("Enter choice (1 or 2): ").strip()
    except:
        choice = "2"  # Default to quick test
    
    if choice == "1":
        main()
    else:
        quick_test()
    
    # Custom test
    print("\nCustom Test:")
    try:
        a = float(input("Enter convection velocity a (e.g., 1.0): ") or "1.0")
        nu = float(input("Enter diffusion coefficient ν (e.g., 0.1): ") or "0.1")
    except:
        a = 1.0
        nu = 0.1
    
    # Compute and plot
    x, u_num = solve_convection_diffusion(a=a, nu=nu)
    u_exact = exact_solution(x, t=1.0, a=a, nu=nu)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, u_num, 'b-', linewidth=2, label='Numerical')
    plt.plot(x, u_exact, 'r--', linewidth=2, alpha=0.7, label='Exact')
    plt.xlabel('x')
    plt.ylabel('u(x,T)')
    plt.title(f'Custom Test: a={a}, ν={nu}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Display error
    error = np.max(np.abs(u_num - u_exact))
    plt.text(0.02, 0.95, f'Maximum absolute error: {error:.2e}', 
            transform=plt.gca().transAxes, fontsize=12,
            verticalalignment='top', 
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.show()