from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

class Su2Grid:
    def __init__( self ):
        self.x = []
        self.y = []
        self.z = []
        self.npoints = -1
        self.nelements = -1
        self.etypes = []
        self.elements = []
        self.etype_map = {}
        self.init_vtkmap()
        self.xmin = -1;
        self.ymin = -1;
        self.zmin = -1;
        self.xmax = 1;
        self.ymax = 1;
        self.zmax = 1;
        self.su2filename = "mesh_RAE2822_turb.su2"
        print( "su2filename = ", self.su2filename )
    def init_vtkmap(self):
        vtk_vertex = 1
        vtk_line = 3
        vtk_triangle = 5
        vtk_quadrilateral = 9
        vtk_tetrahedron = 10
        vtk_hexahedron = 12
        vtk_prism = 13
        vtk_pyramid = 14
        self.etype_map[ vtk_vertex ] = 1
        self.etype_map[ vtk_line ] = 2
        self.etype_map[ vtk_triangle ] = 3
        self.etype_map[ vtk_quadrilateral ] = 4
        self.etype_map[ vtk_tetrahedron ] = 4
        self.etype_map[ vtk_hexahedron ] = 8
        self.etype_map[ vtk_prism ] = 6
        self.etype_map[ vtk_pyramid ] = 5
        print("etype_map=", self.etype_map)
    def ReadPoints( self, f ):#f, npoints, x, y, z ):
        count = 0
        while count < self.npoints :
            count += 1
            line = f.readline()
            xm,ym,zm,id = line.split()
            #print("xm,ym,zm,id=",xm,ym,zm,id)
            self.x.append( float( xm ) )
            self.y.append( float( ym ) )
            self.z.append( float( zm ) )
        self.xmin = min( self.x )
        self.xmax = max( self.x )
        self.ymin = min( self.y )
        self.ymax = max( self.y )
        self.zmin = min( self.z )
        self.zmax = max( self.z )
    def ReadElements( self, f ):#, nelements, elements, etypes ):
        count = 0
        while count < self.nelements :
            count += 1
            line = f.readline()
            print("line=", line)
            str_arr = line.split()
            print("str_arr=", str_arr)
            etype = int( str_arr[0] )
            nvertex = self.etype_map[ etype ]
            print("etype=", etype)
            print("nvertex=", nvertex)
            elem = []
            for index in range( nvertex ):
                print( " index = ", index )
                elem.append( int( str_arr[index+1] ))
            print("elem=", elem)
            self.elements.append( elem )
            self.etypes.append( etype )
            #print("elements=", elements)

    def ReadSu2Grid( self, filename ):
        with open(filename, 'r', encoding='utf-8-sig') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                print( line )
                keyword,value = line.split('=')
                keyword.strip()
                value.strip()
                print( "keyword = ", keyword, " value = ", value )
                if keyword == "NDIME" :
                    dim = value
                elif keyword == "NPOIN" :
                    print( "NPOIN = ", value )
                    self.npoints = int(value)
                    self.ReadPoints( f )
                elif keyword == "NELEM" :
                    print( "NELEM = ", value )
                    self.nelements = int(value)
                    self.ReadElements( f )
                else :
                    print( "otherwise " )
                    break
        f.close()
    def PlotSimple( self ):
        my_quad = Polygon([(0,0), (1,0), (1,1), (0,1)])
        my_tri  = Polygon([(0,1), (1,1), (0.5,2),])

        fig, ax = plt.subplots(1,1)

        a = ax.add_patch( my_tri )
        a.set_fill(False)
        b = ax.add_patch( my_quad )
        b.set_fill(False)


        plt.ylim(-0.05,2.05)
        plt.xlim(-0.55,1.55)

        plt.show()
    def Plot( self ):
        fig, ax = plt.subplots(1,1)
        fig.canvas.set_window_title('OneFLOW CFD Su2Grid Plot Test')
        fig.suptitle('SU2 Rae2822 Grid')
        count = 0        
        while count < self.nelements :
            elem_id = count
            count += 1
            etype = self.etypes[ elem_id ]
            nvertex = self.etype_map[ etype ]
            elem = self.elements[ elem_id ]
            pmlist = []
            for index in range( nvertex ):
                pid = elem[ index ]
                pm = []
                xm = self.x[ pid ]
                ym = self.y[ pid ]
                zm = self.z[ pid ]
                pm.append( xm )
                pm.append( ym )
                pmlist.append( pm )
            #print("pmlist=", pmlist)
            #input()
            my_polygon = Polygon(pmlist)
            var = ax.add_patch( my_polygon )
            var.set_fill(False)

        #plt.ylim(-10.0,10.0)
        #plt.xlim(-10.0,10.0)
        plt.xlim(self.xmin,self.xmax)
        plt.ylim(self.ymin,self.ymax)

        plt.show()        

    def Run( self ):
        self.ReadSu2Grid( self.su2filename )
        self.Plot()
if __name__ == '__main__':
    su2Grid = Su2Grid()
    su2Grid.Run()
