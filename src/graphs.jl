using PyCall, PyPlot

mlines = pyimport("matplotlib.lines")
animation = pyimport("matplotlib.animation")
celluloid = pyimport("celluloid")
dt = pyimport("datetime")
pushfirst!(PyVector(pyimport("sys")."path"), "../src");
machinery = pyimport("importlib.machinery")
loader = machinery.SourceFileLoader("module.name","../src/graphs_py.py")
graph_py = loader.load_module("module.name")