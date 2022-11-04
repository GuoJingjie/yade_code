
O.bodies.append(sphere((0,0,0),0.5))
O.bodies.append(sphere((1,0,0),0.5))
O.bodies.append(sphere((0,1,0),0.5))
O.bodies.append(sphere((0,0,1),0.5))
 
from yade import qt
qt.Controller()
v = qt.View()
rr = qt.Renderer()
graph = GlExtra_AlphaGraph()
rr.extraDrawers = [graph]
graph
