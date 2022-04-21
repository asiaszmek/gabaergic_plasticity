#!/usr/bin/env python
from collections import OrderedDict
import matplotlib.pyplot as plt
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro

def nano_molarity(N, V):
    return 10 * N / V / NA


def pico_sd(N, S):
    return 10 * N / S / NA


def get_grid_list(My_file):
    return np.array(My_file['model']['grid'])


def get_times(My_file, trial='trial0', output="__main__"):
    return np.array(My_file[trial]['output'][output]['times'])

def get_outputs(my_file):
    return my_file['model']['output'].keys()


def get_populations(my_file, trial='trial0', output='__main__'):
    return np.array(my_file[trial]['output'][output]['population'])

def get_all_species(My_file, output="__main__"):
    return [s.decode('utf-8') for s in My_file['model']['output'][output]['species']]

def get_output_regions(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    outputs = {}
    for son in root:
        if son.tag.endswith('OutputScheme'):
            for grandson in son:
                outputs[grandson.get("filename")] = grandson.get("region")
    return outputs

def make_spatial_avg(pop, ind):
    out = zeros((pop.shape[0], len(ind.keys()), pop.shape[2]))
    for i, key in enumerate(ind.keys()):
        first = ind[key][0]
        last = ind[key][-1]
        out[:, i, :] = pop[: ,first:last+1, :].sum(axis=1)
    return out

if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    fig, ax = plt.subplots(1, 1)    
    for fname in sys.argv[1:]:
        my_file = h5py.File(fname, 'r')
        populations = []
        max_len = 0
        for key in my_file:
            if key == "model":
                continue
            populations.append(get_populations(my_file, trial=key))
            if max_len < populations[-1].shape[0]:
                max_len = populations[-1].shape[0]
            
        new_populations = np.zeros((max_len, populations[0].shape[1],
                                    populations[0].shape[2]))
        for i, population in enumerate(populations):
            tot_len = population.shape[0]
            new_populations[:tot_len] += population
        new_populations = new_populations/len(populations)
    
        #output_dict = get_output_regions(my_file)
        grid = get_grid_list(my_file)
        dendrite = []
        dy = 0.2 # um
        for line in grid:
            if line[15] == b"dend":
                dendrite.append(line)
        ind_len = len(dendrite)
        positions = OrderedDict()
        for i in range(ind_len//3):
            positions[i*dy] = [i*3, i*3+1, i*3+2]
                   
        species = get_all_species(my_file)
        
        PP2B_indices = []
        CKp_indices = []
        for i, specie in enumerate(species):
            print(specie)
            if specie == "Epac1cAMP":
               PP2B_indices.append(i)
            elif "CKp" in specie:
                CKp_indices.append(i)
        means_PP2B = []
        for i, population in enumerate(populations):
            mean = np.zeros((population.shape[0]))
            for idx in PP2B_indices:
                mean += population[:,:ind_len, idx].mean(axis=1)
        
            means_PP2B.append(mean)
        means_CKp = []
        for i, population in enumerate(populations):
            mean = np.zeros((population.shape[0]))
            for idx in CKp_indices:
                mean += population[:,:ind_len, idx].mean(axis=1)
            means_CKp.append(mean)

        min_len = min([len(m) for m in means_PP2B])
        new_mean_PP2B = np.zeros((min_len))
        new_mean_CKp = np.zeros((min_len))
        for m in means_PP2B:
            new_mean_PP2B += m[:min_len]/len(means_PP2B)
        for m in means_CKp:
            new_mean_CKp += m[:min_len]/len(means_CKp)

        new_mean = new_mean_PP2B/new_mean_CKp
        ax.plot(new_mean, label=fname)
        ax.legend()
    plt.show()
