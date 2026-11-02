from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

my_quad = Polygon([(0,0), (1,0), (1,1), (0,1)])
my_tri  = Polygon([(0,1), (1,1), (0.5,2),])

fig, ax = plt.subplots(1,1)

a = ax.add_patch( my_tri )
#a.set_fill(False)
b = ax.add_patch( my_quad )
b.set_fill(False)


plt.ylim(-0.05,2.05)
plt.xlim(-0.55,1.55)

plt.show()
