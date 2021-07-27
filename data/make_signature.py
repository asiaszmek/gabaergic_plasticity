#!/usr/bin/env python
from collections import OrderedDict
import matplotlib.pyplot as plt
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro


NA = Avogadro*1e-23
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


def avg_populations(fle):
    populations= []
    max_len = 0
    for key in fle:
        if key == "model":
            continue
        populations.append(get_populations(fle, trial=key))
        if max_len < populations[-1].shape[0]:
            max_len = populations[-1].shape[0]
            
    new_populations = np.zeros((max_len, populations[0].shape[1],
                                populations[0].shape[2]))
    for i, population in enumerate(populations):
        tot_len = population.shape[0]
        new_populations[:tot_len] += population
    new_populations = new_populations/len(populations)
    return new_populations

def get_volumes_positions(fle):
    grid = get_grid_list(my_file)
    dendrite = []
    dy = 0.2 # um
    volumes = {}

    for line in grid:
        if line[15] == b"dend":
            dendrite.append(line)
    ind_len = len(dendrite)
    positions = OrderedDict()
    for i in range(ind_len//3):
        positions[i*dy] = [i*3, i*3+1, i*3+2]
        volumes[i*dy] = (dendrite[i*3][12]+dendrite[i*3+1][12]
                         +dendrite[i*3+2][12])
    return positions, volumes

def get_population_along(specie_idx, new_population, pos, vols, dy=0.2):
    out = np.zeros((new_population.shape[0], len(pos)))
    for i in range(len(pos)):
        p = i*dy
        temp = np.zeros((new_population.shape[0]))
        for idx in pos[p]:
            temp +=  new_population[:, idx, specie_idx]
        out[:, i] =  nano_molarity(temp, vols[p])
    
    return out
                   

if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    fig, ax = plt.subplots(1, 1)
    fig2, ax2 = plt.subplots(1, len(sys.argv)-1)
    pops = []
    for ax_idx, fname in enumerate(sys.argv[1:]):
        my_file = h5py.File(fname, 'r')
        new_populations = avg_populations(my_file)
        positions, volumes = get_volumes_positions(my_file)

        species = get_all_species(my_file)
        ind_len = len(positions)*3    
        PP2B_indices = []
        CKp_indices = []
        for i, specie in enumerate(species):
            if  "PP2BCaMC" in specie:
               PP2B_indices.append(i)
            elif "CKC" in specie or "CKp":
                CKp_indices.append(i)
        means_PP2B = np.zeros((new_populations.shape[0]))

        for idx in PP2B_indices:
            means_PP2B += new_populations[:,:ind_len, idx].mean(axis=1)
        means_CKp = np.zeros((new_populations.shape[0]))
        for idx in CKp_indices:
            means_CKp += new_populations[:,:ind_len, idx].mean(axis=1)

        new_mean = means_CKp/means_PP2B
        ax.plot(new_mean, label=fname)
        ax.legend()
        ca_idx = species.index("Ca")
        ca_dend = get_population_along(ca_idx, new_populations,
                                              positions,
                                              volumes, dy=0.2)
        pops.append(ca_dend)
    min_val = min([ca.min() for ca in pops])
    max_val = max([ca.max() for ca in pops])
    for ix, ax in enumerate(ax2):
        show = ax.imshow(pops[ix], origin="lower",
                         interpolation="none", aspect="auto", vmin=min_val,
                         vmax=max_val, cmap="Reds")
        
    fig2.colorbar(show)
        
    plt.show()
