import matplotlib.pyplot as plt
import numpy as np

class Ghost:
    def __init__(self, xstart, dx, ncells, lr):
        self.lr = lr
        self.xstart = xstart
        self.dx = dx
        self.ncells = ncells
        self.nnodes = self.ncells + 1
    def generate_mesh(self):
        self.x   = np.zeros(self.nnodes, dtype=np.float64)
        self.y   = np.zeros(self.nnodes, dtype=np.float64)
        self.xcc = np.zeros(self.ncells, dtype=np.float64)
        self.ycc = np.zeros(self.ncells, dtype=np.float64)
        for i in range(0, self.nnodes):
            self.x[i] = self.xstart + i * self.dx
            self.y[i] = 0
            
        for i in range(0, self.ncells):
            self.xcc[i] = 0.5*(self.x[i]+self.x[i+1])
            self.ycc[i] = 0
            
    def printinfo(self):
        print(f"Ghost ncells={self.ncells}")
        print(f"Ghost nnodes={self.nnodes}")
        print(f"Ghost xstart={self.xstart}")
        print(f"Ghost lr={self.lr}")
        print(f"Ghost dx={self.dx}")
        print(f"Ghost x={self.x}")
        print(f"Ghost xcc={self.xcc}")
        
        
    def plot(self):
        y0 = 0
        ytext = y0 - 0.1
        for i in range(0, self.ncells): 
            if self.lr == "L":
                mystr = f"${0-i}$"
            else:
                mystr = f"$N+{i+1}$"
            plt.text(self.xcc[i], ytext, mystr, fontsize=12, ha='center')
            
        plt.scatter(self.xcc, self.ycc, s=50, facecolor='red', edgecolor='black', linewidth=1)            
    

class Mesh:
    def __init__(self):
        self.ncells = 9
        self.nnodes = self.ncells + 1
        self.xmin  = 0.0
        self.xmax  = 1.0
        
    def generate_mesh(self):
        self.dx  = (self.xmax-self.xmin) / self.ncells
        print(f"self.dx={self.dx}")
        self.x   = np.zeros(self.nnodes, dtype=np.float64)
        self.y   = np.zeros(self.nnodes, dtype=np.float64)
        self.xcc = np.zeros(self.ncells, dtype=np.float64)
        self.ycc = np.zeros(self.ncells, dtype=np.float64)
        
        for i in range(0, self.nnodes):
            self.x[i] = self.xmin + self.dx*(i)
            self.y[i] = 0
    
        for i in range(0, self.ncells):
            self.xcc[i] = 0.5*(self.x[i]+self.x[i+1])
            self.ycc[i] = 0
            
        #print(f"self.x={self.x}")
        self.nghosts = 2
            
        self.ghost_mesh_left = Ghost(self.x[0],-self.dx,self.nghosts,"L")
        self.ghost_mesh_right = Ghost(self.x[-1],self.dx,self.nghosts,"R")
        
        self.ghost_mesh_left.generate_mesh()
        self.ghost_mesh_right.generate_mesh()
    def printinfo(self):
        print(f"ncells={self.ncells}")
        print(f"nnodes={self.nnodes}")
        print(f"xmin,xmax={self.xmin,self.xmax}")
        self.ghost_mesh_left.printinfo()
        self.ghost_mesh_right.printinfo()
        
    def plot_vertical_interface_lines(self):
        dy = 0.1 * self.dx
        for i in range(0, self.nnodes):
            xm = self.x[i]
            ym = self.y[i]
            coef = 1
            if i == 0 or i== self.nnodes -1:
                coef = 2
            plt.plot([xm, xm], [ym-coef*dy, ym+coef*dy], 'k-')  # 绘制垂直线        
        
    def plot(self):
        plt.scatter(self.xcc, self.ycc, s=50, facecolor='black', edgecolor='black', linewidth=1)
        plt.plot(self.x, self.y, 'k-', linewidth=1)
        dy = 0.1 * self.dx
        self.plot_vertical_interface_lines()
        """
        for i in range(0, self.nnodes):
            xm = self.x[i]
            ym = self.y[i]
            coef = 1
            if i == 0 or i== self.nnodes -1:
                coef = 2
            plt.plot([xm, xm], [ym-coef*dy, ym+coef*dy], 'k-')  # 绘制垂直线
        """
        plt.text(self.x[0], self.y[0]+3*dy, r'$x=a$', fontsize=12, ha='center')
        plt.text(self.x[-1], self.y[0]+3*dy, r'$x=b$', fontsize=12, ha='center')            
        
        self.ghost_mesh_left.plot()
        self.ghost_mesh_right.plot()


def plot_cfd_line( x_points, y0 ):
    # 绘制除中间 5 个点和特定边缘点外的其他点 (内部红色，边缘黑色)
    edge_points_red = np.concatenate([x_points[:2], x_points[:-2]])
    plt.scatter(edge_points_red, np.full_like(edge_points_red, y0), s=100, facecolor='red', edgecolor='black', linewidth=1)
    
    # 绘制左侧第三点 (i=-4) 和右侧第三点 (i=4) 为纯黑色点
    special_black_points = np.array([-4, 4])
    plt.scatter(special_black_points, np.full_like(special_black_points, y0), s=100, facecolor='black', edgecolor='black', linewidth=1)
    
    # 绘制中间 6 个点 (i=-2, -1, 0, 1, 2, 3)
    #middle_points = x_points[3:8]
    middle_points = x_points[3:9]
    plt.scatter(middle_points, np.full_like(middle_points, y0), s=100, facecolor='black', edgecolor='black', linewidth=1)
    
    # 绘制中间 6 个点的黑实线连接
    plt.plot(middle_points, np.full_like(middle_points, y0), 'k-', linewidth=1)
    
    # 添加左起第三点和第四点之间的分段连线（-4到-2）
    plot_mixed_line(-4,-2)
    
    # 添加右起第三点和第四点之间的分段连线（2到4）
    plot_mixed_line(2,4)
    
def plot_cfd_figure():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif=['Times New Roman'])
    
    plt.figure(figsize=(12, 4))

    #x_points = np.array([-20,-10, -7, -6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6], dtype=np.float64)
    inner_points = 7
    bc_points_left = 2
    bc_points_right = 2
    points_max = inner_points + bc_points_right
    x_points = np.arange(-bc_points_left, points_max+1, 1, dtype=np.float64)  # 终止值（开区间），步长1
    print(x_points)  # 输出：[-2. -1.  0.  1.  2.]   
    mesh = Mesh()
    mesh.generate_mesh()
    mesh.printinfo()
    mesh.plot()

    y0 = 0
    

    # Key: Set symmetric axis limits
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1, 1)
    
    plt.axis('equal')
    plt.axis('off')
    
    plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
    plt.show()

if __name__ == '__main__':
    plot_cfd_figure()