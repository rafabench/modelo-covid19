import wquantiles
from matplotlib.cbook import violin_stats
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
import numpy as np
from matplotlib.ticker import EngFormatter, PercentFormatter
import datetime as dt

def preliminar(ax=None, pad=0.3):
    if ax is None:
        ax = plt.gca()
    plt.text(0.5, 0.6, "preliminar",
             alpha=0.2, size=50,
             ha="center", va="top", transform=ax.transAxes,
             bbox=dict(boxstyle="square", pad=pad,
                       ec=(1., 0.8, 0.8, 0.2),
                       fc=(1., 0.8, 0.8, 0.2),
                       )
            )
            
def plot_daily(serie, data, real_data = True, ax = None, label = "Modelo", blur = 1.0, cor_serie = "C0" ,cor_dados = "C1x-", legend = True, cut_day = None):
    n_pts = len(serie)
    start = dt.datetime.strptime("29-03-2020", "%d-%m-%Y")
    then = start + dt.timedelta(days=n_pts-1)
    days = mdates.drange(start,then,dt.timedelta(days=1))
    if ax == None:
        ax = plt.gca()
    ax.plot(days, serie[1:] - serie[:-1], cor_serie, label=label, alpha = blur)
    if cut_day != None:
        new_cut_day = start + dt.timedelta(days=cut_day)
        ax.axvline(x=new_cut_day, color="grey", linestyle="--")
    if real_data:
        ax.plot(days, data, cor_dados)
    ax.grid()
    ax.set_xlabel("Dias")
    ax.yaxis.set_major_formatter(EngFormatter())
    if legend:
        ax.legend(loc = 1)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m"))

def plot_accum(serie, data, real_data = True, ax = None, label = "Modelo", blur = 1.0, cor_serie = "C0", cor_dados = 'C1x-', eng_fmt = True, legend = True, cut_day = None):
    n_pts = len(serie)
    start = dt.datetime.strptime("29-03-2020", "%d-%m-%Y")
    then = start + dt.timedelta(days=n_pts)
    days = mdates.drange(start,then,dt.timedelta(days=1))
    if ax == None:
        ax = plt.gca()
    ax.plot(days, serie, cor_serie, label = label, alpha = blur)
    if cut_day != None:
        new_cut_day = start + dt.timedelta(days=cut_day)
        ax.axvline(x=new_cut_day, color="grey", linestyle="--")
    if real_data:
        ax.plot(days, data, cor_dados, label="Dados")
    ax.grid()
    ax.set_xlabel("Dias")
    if eng_fmt:
        ax.yaxis.set_major_formatter(EngFormatter())
    if legend:
        ax.legend()
    
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m"))
    return

def vdensity_with_weights(weights):
    ''' Outer function allows innder function access to weights. Matplotlib
    needs function to take in data and coords, so this seems like only way
    to 'pass' custom density function a set of weights '''

    def vdensity(data, coords):
        ''' Custom matplotlib weighted violin stats function '''
        # Using weights from closure, get KDE fomr statsmodels
        weighted_cost = sm.nonparametric.KDEUnivariate(data)
        weighted_cost.fit(fft=False, weights=weights)

        # Return y-values for graph of KDE by evaluating on coords
        return weighted_cost.evaluate(coords)
    return vdensity

def custom_violin_stats(data, weights):
    # Get wquantiles median and mean (using wquantiles module for median)
    median = wquantiles.quantile_1D(data, weights, 0.5)
    mean, sumw = np.ma.average(data, weights=list(weights), returned=True)

    # Use matplotlib violin_stats, which expects a function that takes in data and coords
    # which we get from closure above
    results = violin_stats(data, vdensity_with_weights(weights))

    # Update result dictionary with our updated info
    results[0][u"mean"] = mean
    results[0][u"median"] = median

    # No need to do this, since it should be populated from violin_stats
    # results[0][u"min"] =  np.min(data)
    # results[0][u"max"] =  np.max(data)

    return results

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def graph_series(xs,ws,ts, ax = None):
    d = dt.date(2020, 3, 29)
    if ax == None:
        ax = plt.gca()
    count = 1
    for x,w in zip(xs,ws):
        vpstats = custom_violin_stats(x, w)
        vplot = ax.violin(vpstats, [count], vert=True, showmeans=True, showextrema=True,
                      showmedians=True)
        count += 1
    dates = [d + dt.timedelta(days=t) for t in ts]
    labels = [new_d.strftime("%d/%m/%Y") for new_d in dates]
    set_axis_style(ax, labels)
    ax.yaxis.set_major_formatter(EngFormatter())
    return

def graph_rt(xs,ws,ts):
    d = dt.date(2020, 3, 29)
    fig, ax = plt.subplots(figsize=(16,8))
    count = 1
    for x,w in zip(xs,ws):
        vpstats = custom_violin_stats(x, w)
        vplot = ax.violin(vpstats, [count], vert=True, showmeans=True, showextrema=True,
                      showmedians=True)
        count += 1
    dates = [d + dt.timedelta(days=t) for t in ts]
    labels = [new_d.strftime("%d/%m/%Y") for new_d in dates]
    ax.plot(range(0,len(ts)+2),[1 for i in range(0,len(ts)+2)], "k--", alpha=.5)
    set_axis_style(ax, labels)
    ax.set_title("Série da taxa básica de reprodução")
    return

def graph_deaths(xs,ws,ts, day, labels, x_label, reversed = True):
    d = dt.date(2020, 3, 29)
    fig, ax = plt.subplots(figsize=(16,8))
    if reversed:
        count = len(ts)
    else:
        count = 1
    for x,w in zip(xs,ws):
        vpstats = custom_violin_stats(x, w)
        vplot = ax.violin(vpstats, [count], vert=True, showmeans=True, showextrema=True,
                      showmedians=True)
        if reversed:
            count -= 1
        else:
            count += 1
    if reversed:
        labels = labels[::-1]

    set_axis_style(ax, labels)
    #ax.set_yscale("log")
    ax.yaxis.set_major_formatter(EngFormatter())
    new_d = d + dt.timedelta(days=day)
    ax.set_title('Número de mortes até o dia ' + new_d.strftime("%d/%m/%Y"))
    ax.set_xlabel(x_label)
    ax.set_ylabel('Número de mortes')
    #ax.grid()
    plt.subplots_adjust(bottom=0.15, wspace=5)
    return