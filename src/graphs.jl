using PyPlot, PyCall

ticker = PyCall.pyimport("matplotlib.ticker")
mdates = PyCall.pyimport("matplotlib.dates")
dt     = PyCall.pyimport("datetime");

function plot_accum(serie, data;real_data = true, ax = nothing, m_label = "Modelo", blur = 1.0, cor_dados = "x-", eng_fmt = true)
    n_pts = length(serie)
    start = dt.datetime.strptime("29-03-2020", "%d-%m-%Y")
    then = start + dt.timedelta(days=n_pts)
    days = mdates.drange(start,then,dt.timedelta(days=1));
    if ax == nothing
        ax = PyPlot.gca()
    end
    ax.plot(days, serie, label = m_label, alpha = blur)
    if real_data
        ax.plot(days, data, label="Dados", cor_dados)
    end
    ax.grid()
    ax.set_xlabel("Dias")
    if eng_fmt
        ax.yaxis.set_major_formatter(ticker.EngFormatter())
    end
    ax.legend()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m"))
    return
end

function plot_daily(serie, data;real_data = true, ax = nothing, label = "Modelo", blur = 1.0, cor_dados = "x-")
    n_pts = length(serie)
    start = dt.datetime.strptime("29-03-2020", "%d-%m-%Y")
    then = start + dt.timedelta(days=n_pts)
    days = mdates.drange(start,then,dt.timedelta(days=1));
    if ax == nothing
        ax = PyPlot.gca()
    end
    ax.plot(days, serie[2:end] .- serie[1:end-1], label=label, alpha = blur)
    if real_data
        ax.plot(days, data, label="Dados", cor_dados)
    end
    ax.grid()
    ax.set_xlabel("Dias")
    ax.yaxis.set_major_formatter(ticker.EngFormatter())
    ax.legend()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m"))
end

pushfirst!(PyVector(pyimport("sys")."path"), "../src");
machinery = pyimport("importlib.machinery")
loader = machinery.SourceFileLoader("module.name","../src/graphs_violin.py")
graph_violin = loader.load_module("module.name")