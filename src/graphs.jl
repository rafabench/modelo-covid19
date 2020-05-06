using PyCall, PyPlot

pushfirst!(PyVector(pyimport("sys")."path"), "../src");
machinery = pyimport("importlib.machinery")
loader = machinery.SourceFileLoader("module.name","../src/graphs_py.py")
graph_py = loader.load_module("module.name")