'''
Run Sobol sensitivity analysis for the locust model, utilizing Saltelli's sequence.
This uses the SALib library and draws on Python multiprocessing to conduct
simulations in parallel. The resulting data is handled by Pandas.

Author: Christopher Strickland
Date Created: 6/10/2019
'''

import sys, os, time, warnings
from sys import platform
if platform == 'darwin': # OSX backend does not support blitting
    import matplotlib
    matplotlib.use('Qt5Agg')
import pickle
import argparse
from multiprocessing import Pool
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import pandas as pd
from matplotlib import gridspec
import matplotlib.pyplot as plt
import analytic

default_N = os.cpu_count()
parser = argparse.ArgumentParser()
parser.add_argument("-N", type=int, default=1000,
                    help="obtain N*(2D+2) samples from parameter space")
parser.add_argument("-n", "--ncores", type=int,
                    help="number of cores, defaults to {} on this machine".format(default_N))
parser.add_argument("-o", "--filename", type=str, 
                    help="filename to write output to, no extension",
                    default='analysis')

def run_model(N, Rplus, v, loglam, logDelta, eta_alpha, theta_beta, beta,
              gamma, delta):
    '''This is a wrapper around the analytic traveling wave solution based on the
    parameter space that we are examining. It takes in each parameter we are
    testing and calls the solver.
    It then parses the result into whatever we are testing for sensitivity, and
    returns it.
    '''

    # logarithmic lambda, base 10
    lam = 10**loglam
    # logarithmic Delta, base 10
    Delta = 10**logDelta

    # Convert to model parameters in switching dict
    switching = {}
    switching['beta'] = beta
    switching['theta'] = beta*theta_beta
    switching['delta'] = delta

    switching['alpha'] = Delta*beta*switching['theta']/(switching['theta']-beta*eta_alpha)
    switching['eta'] = switching['alpha']*eta_alpha
    switching['gamma'] = gamma

    # Run model
    try:
        speed, R_minus, wave_form, mesh = analytic.main(lam, v, N, Rplus, switching)
    except:
        return (sys.exc_info()[1], None, None, None, None)

    ##### Get peak, width, and skew observables, convert to Rm/Rp, c/v #####
    if not isinstance(speed, Exception):
        ### Peak size
        peak_size = np.max(wave_form)

        ### Width, as std
        mass = np.trapz(wave_form, mesh)
        mean = np.trapz(wave_form*mesh, mesh)/mass

        # integrate to get variance, take sqrt
        width = np.sqrt(np.trapz(wave_form*(mesh-mean)**2, mesh)/mass)

        ### calculate skew
        skew = np.trapz(wave_form*(mesh-mean)**3, mesh)/width**3/mass

        ### Rm/Rp
        RmRp = R_minus/Rplus

        ### cv
        cv = speed/v

        # Traveling wave speed
        # Peak size
        # Width (btwn 12 locusts/m^2)
        # Skew of wave form
        # Resource density behind locusts, asymptotic
        return (cv, peak_size, width, skew, RmRp)

    else:
        # An error was thrown, return as-is for later processing
        return (speed, R_minus, wave_form, mesh)



def main(N, filename, ncores=None, pool=None):
    '''Runs parameter sensitivity on locust model.

    Arguments:
        N: int, used to calculate the number of model realizations
        filename: output file name
        ncores: number of cores we are running on
        pool: multiprocessing pool
    '''

    ### Define the parameter space within the context of a problem dictionary ###
    problem = {
        # number of parameters
        'num_vars' : 10,
        # parameter names
        'names' : ['N', 'Rplus', 'v', 'loglam',
                   'logDelta', 'eta/alpha', 'theta/beta', 'beta',
                   'gamma', 'delta'], 
        # bounds for each corresponding parameter
        'bounds' : [[5000,30000], [120,250], [0.003,0.1], [-10,-4],
                    [-3,0], [0.08,0.8], [1.25,12.5], [0.01,1],
                    [0.0004,0.08], [0.0004,0.08]]
        # # parameter names
        # 'names' : ['N', 'Rplus', 'v', 'loglam',
        #            'beta', 'beta/theta', 'alpha/beta', 'eta/alpha',
        #            'gamma', 'delta'], 
        # # bounds for each corresponding parameter
        # 'bounds' : [[5000,30000], [120,250], [0.0339,0.0532], [-8,-5],
        #             [0.01,1], [0.1,0.8], [1,100], [0.01,0.8],
        #             [0.0004,0.08], [0.0004,0.08]]
    }

    ### Create an N*(2D+2) x num_var matrix of parameter values ###
    param_values = saltelli.sample(problem, N, calc_second_order=True)

    ### Run model ###
    print('Examining the parameter space.')
    if pool is not None:
        if ncores is None:
            poolsize = os.cpu_count()
        else:
            poolsize = ncores
        chunksize = param_values.shape[0]//poolsize
        output = pool.starmap(run_model, param_values, chunksize=chunksize)
    else:
        # This is a bad idea. Say so.
        print('Warning!!! Running in serial only! Multiprocessing pool not utilized.')
        output = []
        for params in param_values:
            output.append(run_model(*params))

    ### Parse and save the output ###
    print('Saving and checking the results for errors...')
    # write data to temporary location in case anything obnoxious happens
    with open("raw_result_data.pickle", "wb") as f:
        result = {'output':output, 'param_values':param_values}
        pickle.dump(result, f)
    # Look for errors
    error_num = 0
    error_places = []
    for n, result in enumerate(output):
        if isinstance(result[0], Exception):
            error_num += 1
            error_places.append(n)
    if error_num > 0:
        print("Errors discovered in output.")
        print("Parameter locations: {}".format(error_places))
        print("Pickling errors...")
        err_output = []
        err_params = []
        for idx in error_places:
            err_output.append(output.pop(idx)) # pull out err output
            err_params.append(param_values[idx,:]) # reference err param values
        with open("err_results.pickle", "wb") as f:
            err_result = {'err_output':err_output, 'err_param_values':param_values}
            pickle.dump(err_result, f)
        print("Saving all other data in HDF5...")
        output = np.array(output)
        # first remove err param values
        param_values = param_values[[ii for ii in range(len(param_values)) if ii not in error_places],:]
        # convert to dataframe
        param_values = pd.DataFrame(param_values, columns=problem['names'])
        assert output.shape[0] == param_values.shape[0]
        store = pd.HDFStore('nonerr_results.h5')
        store['param_values'] = param_values
        store['raw_output'] = pd.DataFrame(output, 
                            columns=['c/v', 'peak', 'width', 'skew', 'Rm/Rp'])
        store.close()
        os.remove('raw_result_data.pickle')
        print("Please review output dump.")
        return
    # Save results in HDF5 as dataframe
    print('Parsing the results...')
    output = np.array(output)
    param_values = pd.DataFrame(param_values, columns=problem['names'])
    # Resave as dataframe in hdf5
    store = pd.HDFStore(filename+'.h5')
    store['param_values'] = param_values
    store['raw_output'] = pd.DataFrame(output, 
                          columns=['c/v', 'peak', 'width', 'skew', 'Rm/Rp'])
    os.remove('raw_result_data.pickle')

    ### Analyze the results and view using Pandas ###
    # Conduct the sobol analysis and pop out the S2 results to a dict
    S2 = {}
    cv_sens = sobol.analyze(problem, output[:,0], calc_second_order=True)
    S2['cv'] = pd.DataFrame(cv_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['cv_conf'] = pd.DataFrame(cv_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    peak_sens = sobol.analyze(problem, output[:,1], calc_second_order=True)
    S2['peak'] = pd.DataFrame(peak_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['peak_conf'] = pd.DataFrame(peak_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    width_sens = sobol.analyze(problem, output[:,2], calc_second_order=True)
    S2['width'] = pd.DataFrame(width_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['width_conf'] = pd.DataFrame(width_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    skew_sens = sobol.analyze(problem, output[:,3], calc_second_order=True)
    S2['skew'] = pd.DataFrame(skew_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['skew_conf'] = pd.DataFrame(skew_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    RmRp_sens = sobol.analyze(problem, output[:,4], calc_second_order=True)
    S2['RmRp'] = pd.DataFrame(RmRp_sens.pop('S2'), index=problem['names'],
                           columns=problem['names'])
    S2['RmRp_conf'] = pd.DataFrame(RmRp_sens.pop('S2_conf'), index=problem['names'],
                           columns=problem['names'])
    # Convert the rest to a pandas dataframe
    cv_sens = pd.DataFrame(cv_sens,index=problem['names'])
    peak_sens = pd.DataFrame(peak_sens,index=problem['names'])
    width_sens = pd.DataFrame(width_sens,index=problem['names'])
    skew_sens = pd.DataFrame(skew_sens,index=problem['names'])
    RmRp_sens = pd.DataFrame(RmRp_sens,index=problem['names'])

    ### Save the analysis ###
    print('Saving...')
    store['cv_sens'] = cv_sens
    store['peak_sens'] = peak_sens
    store['width_sens'] = width_sens
    store['skew_sens'] = skew_sens
    store['RmRp_sens'] = RmRp_sens
    for key in S2.keys():
        store['S2/'+key] = S2[key]
    # Save the bounds
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        store['bounds'] = pd.Series(problem['bounds'])
    store.close()

    # Plot and save
    plot_S1_ST_tbl(cv_sens, peak_sens, width_sens, skew_sens, RmRp_sens,
                   problem['bounds'], False)



def load_data(filename):
    '''Load analysis data from previous run and return for examination'''

    return pd.HDFStore(filename)



def print_max_conf(store):
    '''Print off the max confidence interval for each variable in the store,
    for both first-order and total-order indices'''
    for var in ['cv_sens', 'peak_sens', 'width_sens', 'skew_sens', 'RmRp_sens']:
        print('----------- '+var+' -----------')
        print('S1_conf_max: {}'.format(store[var]['S1_conf'].max()))
        print('ST_conf_max: {}'.format(store[var]['ST_conf'].max()))
        print(' ')



def find_feasible_params(store):
    '''This will return a dataframe with parameters and corresponding 
    observable values that fall within the desired ranges (hard coded)'''
    data = pd.concat([store['param_values'], store['raw_output']], axis=1)
    data1 = data[data['skew']>=1]
    data1 = data1[data1['skew']<=2]
    if data1.shape[0] == 0:
        print('No parameters match skew requirement alone.')
        return data
    data2 = data1[data1['peak']<=4280]
    data2 = data2[data2['peak']>=950]
    if data2.shape[0] == 0:
        print('No parameters match both skew and peak requirements.')
        return data1
    data3 = data2[data2['c/v']*data2['v']<=0.00852]
    data3 = data3[data3['c/v']*data3['v']>=0.000508]
    if data3.shape[0] == 0:
        print('No parameters match skew, peak, and vel requirements.')
        return data2

    return data3


def plot_S1_ST_tbl_from_store(store, show=True, ext='pdf'):
    cv_sens = store['cv_sens']
    peak_sens = store['peak_sens']
    width_sens = store['width_sens']
    skew_sens = store['skew_sens']
    RmRp_sens = store['RmRp_sens']
    bounds = list(store['bounds'])
    
    plot_S1_ST_tbl(peak_sens, width_sens, skew_sens, cv_sens, RmRp_sens,
                   bounds, show, ext)



def plot_S1_ST_tbl(peak_sens, width_sens, skew_sens, cv_sens, RmRp_sens,
                   bounds=None, show=True, ext='pdf', startclr=None):

    # Gather the S1 and ST results
    all_names = ['max density', 'std width', 'skewness', 'c/v', r'$R^-/R^+$']
    all_results = [peak_sens, width_sens, skew_sens, cv_sens, RmRp_sens]
    names = []
    results = []
    for n, result in enumerate(all_results):
        # Keep only the ones actually passed
        if result is not None:
            results.append(result)
            names.append(all_names[n])
    S1 = pd.concat([result['S1'] for result in results[::-1]], keys=names[::-1], axis=1) #produces copy
    ST = pd.concat([result['ST'] for result in results[::-1]], keys=names[::-1], axis=1)

    # Gather the S1 and ST results without skew
    # S1 = pd.concat([cv_sens['S1'], peak_sens['S1'], width_sens['S1'], 
    #                RmRp_sens['S1']], 
    #                keys=['$c$', 'peak', 'std', '$R^-$'], axis=1) #produces copy
    # ST = pd.concat([cv_sens['ST'], peak_sens['ST'], width_sens['ST'], 
    #                RmRp_sens['ST']], 
    #                keys=['$c$', 'peak', 'std', '$R^-$'], axis=1)
    ##### Reorder (manual) #####
    order = ['N', 'Rplus', 'v', 'beta', 'eta/alpha', 'theta/beta' ,'logDelta', 
             'gamma', 'delta', 'loglam']
    bndorder = [0, 1, 2, -3, 5, 6, 4, -2, -1, 3]
    S1 = S1.reindex(order)
    ST = ST.reindex(order)
    if bounds is not None:
        new_bounds = [bounds[ii] for ii in bndorder]
        bounds = new_bounds

    ###### Change to greek, LaTeX #####
    for id in S1.index:
        # handle greek ratios
        if '/' in id:
            togreek = "$\\" + id[:id.find('/')+1] + "\\" + id[id.find('/')+1:] + r"$"
            S1.rename(index={id: togreek}, inplace=True)
        elif id == 'loglam':
            S1.rename(index={id: r'$\log(\lambda)$'}, inplace=True)
        elif id == 'logDelta':
            S1.rename(index={id: r'$\log(\Delta)$'}, inplace=True)
        elif id == 'Rplus':
            S1.rename(index={id: r'$R^{+}$'}, inplace=True)
        # all others
        elif id not in ['N', 'v']:
            S1.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
        
    for id in ST.index:
        if '/' in id:
            togreek = "$\\" + id[:id.find('/')+1] + "\\" + id[id.find('/')+1:] + r"$"
            ST.rename(index={id: togreek}, inplace=True)
        elif id == 'loglam':
            ST.rename(index={id: r'$\log(\lambda)$'}, inplace=True)
        elif id == 'logDelta':
            ST.rename(index={id: r'$\log(\Delta)$'}, inplace=True)
        elif id == 'Rplus':
            ST.rename(index={id: r'$R^{+}$'}, inplace=True)
        elif id not in ['N', 'v']:
            ST.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
    
    ###### Plot ######
    if bounds is not None:
        # setup for table
        fig = plt.figure(figsize=(15, 6))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,.35], wspace=.15, left=0.04,
                               right=0.975, bottom=0.15, top=0.915)
        axes = []
        for ii in range(2):
            axes.append(plt.subplot(gs[ii]))
    else:
        # setup without table
        fig, axes = plt.subplots(ncols=2, figsize=(13, 6))
    # Switch the last two colors so A is red
    # prop_cycle = plt.rcParams['axes.prop_cycle']
    # colors = prop_cycle.by_key()['color']
    # colors = colors[:4]
    # clr = colors[-1]; colors[-1] = colors[-2]; colors[-2] = clr
    
    if startclr is not None:
        # Start at a different color
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        colors = colors[startclr:]
        s1bars = S1.plot.bar(stacked=True, ax=axes[0], rot=0, width=0.8, color=colors)
        s2bars = ST.plot.bar(stacked=True, ax=axes[1], rot=0, width=0.8, color=colors, legend=False)
    else:
        s1bars = S1.plot.bar(stacked=True, ax=axes[0], rot=0, width=0.8)
        s2bars = ST.plot.bar(stacked=True, ax=axes[1], rot=0, width=0.8, legend=False)
    for ax in axes:
        ax.tick_params(axis='x', labelsize=18, rotation=-40) #-25
        ax.tick_params(axis='y', labelsize=14)
        #ax.get_yaxis().set_visible(False)
        ax.set_ylim(bottom=0)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    handles, labels = s1bars.get_legend_handles_labels()
    s1bars.legend(reversed(handles), reversed(labels), loc='upper left', fontsize=16)
    # Create table
    if bounds is not None:
        columns = ('Value Range',)
        rows = list(S1.index)
        # turn bounds into strings of ranges
        cell_text = []
        for bnd in bounds:
            low = str(bnd[0])
            high = str(bnd[1])
            # concatenate, remove leading zeros
            if low != "0" and low != "0.0":
                low = low.lstrip("0")
            if high != "0" and high != "0.0":
                high = high.lstrip("0")
            # raise any minus signs
            if '-' in low:
                low = "\u00AF"+low[1:]
            if '-' in high:
                high = "\u00AF"+high[1:]
            cell_text.append([low+"-"+high])
        tbl_ax = plt.subplot(gs[2])
        the_table = tbl_ax.table(cellText=cell_text, rowLabels=rows, colLabels=columns,
                    loc='center')
        the_table.set_fontsize(18)
        the_table.scale(1,2.3)
        the_table.auto_set_column_width(0)
        tbl_ax.axis('off')
    #plt.tight_layout()
    # reposition table
    pos = tbl_ax.get_position()
    newpos = [pos.x0 + 0.02, pos.y0, pos.width, pos.height]
    tbl_ax.set_position(newpos)
    if show:
        plt.show()
    else:
        fig.savefig("param_sens_{}.{}".format(time.strftime("%m_%d_%H%M"), ext))
    return (fig, axes)



def plot_shape_rel_from_store(store, show=True, ext='pdf'):
    cv_sens = store['cv_sens']
    peak_sens = store['peak_sens']
    width_sens = store['width_sens']
    skew_sens = store['skew_sens']
    RmRp_sens = store['RmRp_sens']
    bounds = list(store['bounds'])
    
    plot_shape_rel(peak_sens, width_sens, skew_sens, cv_sens, RmRp_sens,
                   bounds, show, ext)



def plot_shape_rel(peak_sens, width_sens, skew_sens, cv_sens, RmRp_sens,
                   bounds=None, show=True, ext='pdf'):

    # Gather the shape and rel results into two groups
    shape_results = [peak_sens, width_sens, skew_sens]
    shape_names = ['max density', 'std width', 'skewness']
    rel_names = ['c/v', r'$R^-/R^+$']
    rel_results = [cv_sens, RmRp_sens]

    S1_shape = pd.concat([result['S1'] for result in shape_results[::-1]], keys=shape_names[::-1], axis=1) #produces copy
    ST_shape = pd.concat([result['ST'] for result in shape_results[::-1]], keys=shape_names[::-1], axis=1)
    S1_rel = pd.concat([result['S1'] for result in rel_results[::-1]], keys=rel_names[::-1], axis=1)
    ST_rel = pd.concat([result['ST'] for result in rel_results[::-1]], keys=rel_names[::-1], axis=1)

    ##### Reorder parameters (manual) #####
    order = ['N', 'Rplus', 'v', 'beta', 'eta/alpha', 'theta/beta' ,'logDelta', 
             'gamma', 'delta', 'loglam']
    bndorder = [0, 1, 2, -3, 5, 6, 4, -2, -1, 3]
    S1_shape = S1_shape.reindex(order)
    ST_shape = ST_shape.reindex(order)
    S1_rel = S1_rel.reindex(order)
    ST_rel = ST_rel.reindex(order)
    if bounds is not None:
        new_bounds = [bounds[ii] for ii in bndorder]
        bounds = new_bounds

    ###### Change to greek, LaTeX #####
    for id in S1_shape.index:
        # handle greek ratios
        if '/' in id:
            togreek = "$\\" + id[:id.find('/')+1] + "\\" + id[id.find('/')+1:] + r"$"
            S1_shape.rename(index={id: togreek}, inplace=True)
            ST_shape.rename(index={id: togreek}, inplace=True)
            S1_rel.rename(index={id: togreek}, inplace=True)
            ST_rel.rename(index={id: togreek}, inplace=True)
        elif id == 'loglam':
            S1_shape.rename(index={id: r'$\log(\lambda)$'}, inplace=True)
            ST_shape.rename(index={id: r'$\log(\lambda)$'}, inplace=True)
            S1_rel.rename(index={id: r'$\log(\lambda)$'}, inplace=True)
            ST_rel.rename(index={id: r'$\log(\lambda)$'}, inplace=True)
        elif id == 'logDelta':
            S1_shape.rename(index={id: r'$\log(\Delta)$'}, inplace=True)
            ST_shape.rename(index={id: r'$\log(\Delta)$'}, inplace=True)
            S1_rel.rename(index={id: r'$\log(\Delta)$'}, inplace=True)
            ST_rel.rename(index={id: r'$\log(\Delta)$'}, inplace=True)
        elif id == 'Rplus':
            S1_shape.rename(index={id: r'$R^{+}$'}, inplace=True)
            ST_shape.rename(index={id: r'$R^{+}$'}, inplace=True)
            S1_rel.rename(index={id: r'$R^{+}$'}, inplace=True)
            ST_rel.rename(index={id: r'$R^{+}$'}, inplace=True)
        # all others
        elif id not in ['N', 'v']:
            S1_shape.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
            ST_shape.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
            S1_rel.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
            ST_rel.rename(index={id: r'$\{}$'.format(id)}, inplace=True)
            
    ###### Plot ######
    if bounds is not None:
        # setup for table
        fig = plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(2, 3, figure=fig, width_ratios=[1,1,.35], 
                    wspace=.15, left=0.04, right=0.975, bottom=0.08, top=0.95)
        axes = []
        for ii in range(2):
            for jj in range(2):
                axes.append(plt.subplot(gs[ii,jj]))
    else:
        # setup without table
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(13, 12))

    s1bars_shape = S1_shape.plot.bar(stacked=True, ax=axes[0], rot=0, width=0.8)
    s2bars_shape = ST_shape.plot.bar(stacked=True, ax=axes[1], rot=0, width=0.8, legend=False)
    # Start relative bars at third color
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    colors = colors[3:]
    s1bars_rel = S1_rel.plot.bar(stacked=True, ax=axes[2], rot=0, width=0.8, color=colors)
    s2bars_rel = ST_rel.plot.bar(stacked=True, ax=axes[3], rot=0, width=0.8, color=colors, legend=False)
        

    for ax in axes:
        ax.tick_params(axis='x', labelsize=18, rotation=-40) #-25
        ax.tick_params(axis='y', labelsize=14)
        #ax.get_yaxis().set_visible(False)
        ax.set_ylim(bottom=0)
    axes[0].set_title('First-order indices', fontsize=26)
    axes[1].set_title('Total-order indices', fontsize=26)
    
    # set legends
    for s1bars in [s1bars_shape, s1bars_rel]:
        handles, labels = s1bars.get_legend_handles_labels()
        s1bars.legend(reversed(handles), reversed(labels), loc='upper left', fontsize=16)
    # Create table
    if bounds is not None:
        columns = ('Value Range',)
        rows = list(S1_shape.index)
        # turn bounds into strings of ranges
        cell_text = []
        for bnd in bounds:
            low = str(bnd[0])
            high = str(bnd[1])
            # concatenate, remove leading zeros
            if low != "0" and low != "0.0":
                low = low.lstrip("0")
            if high != "0" and high != "0.0":
                high = high.lstrip("0")
            # raise any minus signs
            if '-' in low:
                low = "\u00AF"+low[1:]
            if '-' in high:
                high = "\u00AF"+high[1:]
            cell_text.append([low+"-"+high])
        tbl_ax = plt.subplot(gs[:,2])
        the_table = tbl_ax.table(cellText=cell_text, rowLabels=rows, colLabels=columns,
                    loc='center')
        the_table.set_fontsize(18)
        the_table.scale(1,2.3)
        the_table.auto_set_column_width(0)
        tbl_ax.axis('off')
    #plt.tight_layout()
    # reposition table
    pos = tbl_ax.get_position()
    newpos = [pos.x0 + 0.02, pos.y0, pos.width, pos.height]
    tbl_ax.set_position(newpos)
    if show:
        plt.show()
    else:
        fig.savefig("param_sens_{}.{}".format(time.strftime("%m_%d_%H%M"), ext))
    return (fig, axes)



def plot_RmRp_from_store(store, peak_trunc=None, down_samp_frac=None, xrange=None,
                         show=True, save_png=False, save_pdf=False, save_tif=False,
                         plot_point=True):
    '''Function to plot the R_minus/R_plus vs. lambda plot'''

    plt.close('all')
    title_size = 28 #18
    label_size = 24 #14, 21
    tick_size = 24 #14
    fig = plt.figure(figsize=(10,7.5)) #8,6
    gs = gridspec.GridSpec(1, 1, left=0.04, right=0.975, bottom=0.15, top=0.915)
    x = store['param_values']['loglam']
    y = store['raw_output']['Rm/Rp']
    c = store['param_values']['logDelta']
    if down_samp_frac is not None:
        sample_len = len(store['param_values']['loglam'])
        down_num = int(down_samp_frac*sample_len)
        idx = np.random.randint(0, len(store['param_values']['loglam']), down_num)
        x = x[idx]
        y = y[idx]
        c = c[idx]
    if peak_trunc is not None:
        plt.scatter(x[store['raw_output']['peak']<=peak_trunc],
                y[store['raw_output']['peak']<=peak_trunc],
                s=1, c=c[store['raw_output']['peak']<=peak_trunc],
                cmap='viridis_r')
        c_ticks = np.linspace(c[store['raw_output']['peak']<=peak_trunc].min(),
                              c[store['raw_output']['peak']<=peak_trunc].max(), 5,
                              endpoint=True)
        plt.title('Fraction of resources remaining,\ntruncated at peak={:,.0f}'.format(peak_trunc),
            fontsize=title_size)
    else:
        plt.scatter(x, y, s=1, c=c, cmap='viridis_r')
        c_ticks = np.linspace(c.min(), c.max(), 5, endpoint=True)
        plt.title('Fraction of resources remaining', fontsize=title_size)
    if plot_point:
        plt.plot(-5,10**-6,'ro')
    clb = plt.colorbar(ticks=c_ticks)
    clb.ax.set_title(r'$\log_{10}(\Delta)$', fontsize=label_size, pad=15)
    clb.ax.tick_params(labelsize=tick_size)
    plt.xlabel(r"$\log_{10}(\lambda)$", fontsize=label_size)
    plt.ylabel(r"$\frac{R^-}{R^+}$", fontsize=label_size+8, rotation=0)
    ax = plt.gca()
    ax.yaxis.set_label_coords(-.1,0.45)
    if xrange is None:
        plt.xlim(x.min(),x.max())
    else:
        plt.xlim(xrange)
    plt.ylim(0,1)
    plt.tick_params(labelsize=tick_size)
    plt.tight_layout()
    if save_png:
        fig.savefig("resource_plot.png", dpi=300)
    if save_pdf:
        fig.savefig("resource_plot.pdf")
    if save_tif:
        fig.savefig("resource_plot.tif", dpi=300)
    if show:
        plt.show()



def plot_cv_from_store(store, peak_trunc=None, down_samp_frac=None, xrange=None,
                       show=True, save_png=False, save_pdf=False, save_tif=False,
                       plot_point=True):
    '''Function to plot the theta/beta vs. c/v plot'''

    plt.close('all')
    title_size = 32 #18
    label_size = 24 #14, 21
    tick_size = 24 #14
    fig = plt.figure(figsize=(10,7.5)) #8,6
    gs = gridspec.GridSpec(1, 1, left=0.04, right=0.975, bottom=0.15, top=0.915)
    x = store['param_values']['loglam']
    y = store['raw_output']['c/v']
    c = store['param_values']['logDelta']
    if down_samp_frac is not None:
        sample_len = len(store['param_values']['loglam'])
        down_num = int(down_samp_frac*sample_len)
        idx = np.random.randint(0, len(store['param_values']['loglam']), down_num)
        x = x[idx]
        y = y[idx]
        c = c[idx]
    if peak_trunc is not None:
        plt.scatter(x[store['raw_output']['peak']<=peak_trunc],
                y[store['raw_output']['peak']<=peak_trunc],
                s=1, c=c[store['raw_output']['peak']<=peak_trunc],
                cmap='viridis_r')
        c_ticks = np.linspace(c[store['raw_output']['peak']<=peak_trunc].min(),
                              c[store['raw_output']['peak']<=peak_trunc].max(), 5,
                              endpoint=True)
        plt.title('Fraction of collective speed,\ntruncated at peak={:,.0f}'.format(peak_trunc),
            fontsize=title_size)
    else:
        plt.scatter(x, y, s=1, c=c, cmap='viridis_r')
        c_ticks = np.linspace(c.min(), c.max(), 5, endpoint=True)
        plt.title('Fraction of collective speed', fontsize=title_size)
    if plot_point:
        plt.plot(-5,0.1325,'ro')
    clb = plt.colorbar(ticks=c_ticks)
    clb.ax.set_title(r'$\log_{10}(\Delta)$', fontsize=label_size, pad=15)
    clb.ax.tick_params(labelsize=tick_size)
    plt.xlabel(r"$\log_{10}(\lambda)$", fontsize=label_size)
    plt.ylabel(r"$\frac{c}{v}$", fontsize=label_size+10, rotation=0)
    ax = plt.gca()
    ax.yaxis.set_label_coords(-.1,0.47)
    if xrange is None:
        plt.xlim(x.min(),x.max())
    else:
        plt.xlim(xrange)
    #plt.xscale("log")
    plt.ylim(0,y.max())
    plt.tick_params(labelsize=tick_size)
    plt.tight_layout()
    if save_png:
        fig.savefig("speed_plot.png", dpi=300)
    if save_pdf:
        fig.savefig("speed_plot.pdf")
    if save_tif:
        fig.savefig("speed_plot.tif", dpi=300)
    if show:
        plt.show()



def plot_skew_from_store(store, peak_trunc=1e4, down_samp_frac=None, xrange=None,
                         show=True, save_png=False, save_pdf=False, plot_point=True):
    '''Function to plot the R_minus/R_plus vs. lambda plot'''

    plt.close('all')
    # title_size = 18
    # label_size = 14
    # tick_size = label_size
    title_size = 36
    label_size = 24
    tick_size = 24
    fig = plt.figure(figsize=(10,7.5))
    gs = gridspec.GridSpec(1, 1, left=0.04, right=0.975, bottom=0.15, top=0.915)
    x = store['param_values']['loglam']
    y = store['raw_output']['skew']
    c = store['param_values']['logDelta']
    if down_samp_frac is not None:
        sample_len = len(store['param_values']['loglam'])
        down_num = int(down_samp_frac*sample_len)
        idx = np.random.randint(0, len(store['param_values']['loglam']), down_num)
        x = x[idx]
        y = y[idx]
        c = c[idx]
    if peak_trunc is not None:
        plt.scatter(x[store['raw_output']['peak']<=peak_trunc],
                y[store['raw_output']['peak']<=peak_trunc],
                s=1, c=c[store['raw_output']['peak']<=peak_trunc],
                cmap='viridis_r')
        c_ticks = np.linspace(c[store['raw_output']['peak']<=peak_trunc].min(),
                              c[store['raw_output']['peak']<=peak_trunc].max(), 5,
                              endpoint=True)
        plt.title('max density '+r'$\leq$'+'{:,.0f}   '.format(peak_trunc), fontsize=title_size)
    else:
        plt.scatter(x, y, s=1, c=c, cmap='viridis_r')
        c_ticks = np.linspace(c.min(), c.max(), 5, endpoint=True)
        plt.title('skewness plot', fontsize=title_size)
    if plot_point:
        plt.plot(-5,1.78,'ro')
    clb = plt.colorbar(ticks=c_ticks)
    clb.ax.set_title(r"$\log_{10}(\Delta)$", fontsize=label_size, pad=15)
    clb.ax.tick_params(labelsize=tick_size)
    plt.xlabel(r'$\log_{10}(\lambda)$', fontsize=label_size)
    plt.ylabel("skewness", fontsize=label_size, rotation=90)
    ax = plt.gca()
    ax.yaxis.set_label_coords(-0.05,0.56)
    if xrange is None:
        plt.xlim(x.min(),x.max())
    else:
        plt.xlim(xrange)
    #plt.xscale("log")
    #plt.ylim(-8,-5)
    plt.tick_params(labelsize=tick_size)
    plt.tight_layout()
    if save_png:
        fig.savefig("skew_plot.png", dpi=300)
    if save_pdf:
        fig.savefig("skew_plot.pdf")
    if show:
        plt.show()



if __name__ == "__main__":
    args = parser.parse_args()
    if args.ncores is None:
        with Pool() as pool:
            main(args.N, args.filename, args.ncores, pool)
    elif args.ncores > 1:
        with Pool(args.ncores) as pool:
            main(args.N, args.filename, args.ncores, pool)
    else:
        main(args.N, args.filename)
