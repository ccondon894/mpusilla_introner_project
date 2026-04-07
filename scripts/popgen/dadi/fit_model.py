# Pickle is used to save variables as files for future use
import pickle
# NLopt is the optimization libary dadi uses
import nlopt
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
import dadi

datafile = '../data_prep/mpusilla.snps.4d.notMT.nomissing.filtered.vcf.gz' 
dd = dadi.Misc.make_data_dict_vcf(datafile, 'popfile.txt')
pop_ids = ["intronerful","intronerless"]
ns = [11,2]

data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=False)
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

# demographic model
demo_model = dadi.Demographics2D.split_mig

# Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

# Define starting parameters
params = [0.7, 0.7, 12, 0.005]
lower_bounds = [1e-2, 1e-2, 1e-3, 1e-4]
upper_bounds = [3, 3, 30, 1]

def fit_and_print ( fs ) :
    ## fit our model
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                                  lower_bound=lower_bounds)
    popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l,
                                        lower_bound=lower_bounds,
                                        upper_bound=upper_bounds,
                                        algorithm=nlopt.LN_BOBYQA,
                                        maxeval=5000, verbose=10)

    ##
    # Calculate the synonymous theta
    model_fs = demo_model_ex(popt, ns, pts_l)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)

    # Write results to fid
    res = [ll_model] + list(popt) + [theta0]
    try:
      fid = open('model_fits.4d.txt','a')
    except:
      fid = open('model_fits.4d.txt','w')
    fid.write('\t'.join([str(ele) for ele in res])+'\n')
    fid.close()

## main model
fit_and_print( data_fs ) 

## now bootstrap it
chunks = dadi.Misc.fragment_data_dict(dd, 250000)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, 100, pop_ids, ns, polarized=False)
for b in boots :
   fit_and_print( b )