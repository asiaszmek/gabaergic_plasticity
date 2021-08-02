import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
keys = ["10 uM Glu", "50 uM NMDA", "100 uM Glu"] 
filenames_Fluo = {
    "100 uM Glu":
    [
        "model_no_ER_250_nM_spine_Glu_100_uM_Fluo3_trial0_total.txt",
        "model_no_ER_250_nM_spine_Glu_100_uM_Fluo3_trial1_total.txt",
        "model_no_ER_250_nM_spine_Glu_100_uM_Fluo3_trial2_total.txt",
        "model_no_ER_250_nM_spine_Glu_100_uM_Fluo3_trial3_total.txt",
    ],
    "10 uM Glu":
    [
        "model_no_ER_150_nM_spine_Glu_10_uM_Fluo3_trial0_total.txt",
        "model_no_ER_150_nM_spine_Glu_10_uM_Fluo3_trial1_total.txt",
        "model_no_ER_150_nM_spine_Glu_10_uM_Fluo3_trial2_total.txt",
        "model_no_ER_150_nM_spine_Glu_10_uM_Fluo3_trial3_total.txt",
    ],
    "50 uM NMDA":
    [
        "model_no_ER_150_nM_spine_NMDA_Fluo3_trial0_total.txt",
        "model_no_ER_150_nM_spine_NMDA_Fluo3_trial1_total.txt",
        "model_no_ER_150_nM_spine_NMDA_Fluo3_trial2_total.txt",
        "model_no_ER_150_nM_spine_NMDA_Fluo3_trial3_total.txt",
    ]
}
colors = {
    "100 uM Glu": "tab:orange",
    "50 uM NMDA": "tab:blue",
    "10 uM Glu": "tab:green",
}

filenames_Ca = {
    "100 uM Glu":
    [
        "model_no_ER_250_nM_spine_Glu_100_uM_trial0_total.txt",
        "model_no_ER_250_nM_spine_Glu_100_uM_trial1_total.txt",
        "model_no_ER_250_nM_spine_Glu_100_uM_trial2_total.txt",
        "model_no_ER_250_nM_spine_Glu_100_uM_trial3_total.txt",
    ],
    "10 uM Glu":
    [
        "model_no_ER_150_nM_spine_Glu_10_uM_trial0_total.txt",
        "model_no_ER_150_nM_spine_Glu_10_uM_trial1_total.txt",
        "model_no_ER_150_nM_spine_Glu_10_uM_trial2_total.txt",
        "model_no_ER_150_nM_spine_Glu_10_uM_trial3_total.txt",
    ],
    "50 uM NMDA":
    [
        "model_no_ER_150_nM_spine_NMDA_trial0_total.txt",
        "model_no_ER_150_nM_spine_NMDA_trial1_total.txt",
        "model_no_ER_150_nM_spine_NMDA_trial2_total.txt",
        "model_no_ER_150_nM_spine_NMDA_trial3_total.txt",
    ]
}

def get_mean(file_list, specie):
    data = []

    for fname in file_list:
        f = open(fname, "r")
        header = f.readline().split()
        datas = np.loadtxt(f)
        time = datas[:, 0]
        idx = header.index(specie)
        data.append(datas[:, idx])

    min_len = min([len(d) for d in data])
    out = np.zeros((min_len, len(data)))
    for i, d in enumerate(data):
        out[:, i] = np.array(d[:min_len])
    mean_signal = out.mean(axis=1)
    return mean_signal, time[0:min_len]


if __name__ == "__main__":
    t_stim = 100000
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    for key in keys:
        mean_fluo, time = get_mean(filenames_Fluo[key], "Fluo3Ca")
        dt = time[1] - time[0]
        fluo_basal = mean_fluo[int((t_stim-60000)//dt):int(t_stim//dt)].mean()
        fluo_signal = (mean_fluo-fluo_basal)/fluo_basal

        time = (time - t_stim)/1000/60 # convert to minutes
        ax[0].plot(time, fluo_signal, color=colors[key], marker="o", label=key,
                linestyle = 'None', markerfacecolor='None',
                markeredgecolor=colors[key])
        mean_ca, time = get_mean(filenames_Ca[key], "Ca")
        dt = time[1] - time[0]
        time = (time - t_stim)/1000/60 # convert to minutes
        ax[1].plot(time, mean_ca, color=colors[key], marker="d", label=key,
                 linestyle = 'None', markerfacecolor='None',
                 markeredgecolor=colors[key])
    ax[0].set_ylabel("Fluorescence", fontsize=14)
    ax[1].set_ylabel("Ca concentration (nM)", fontsize=14)
    ax[0].set_xlabel("time (min)", fontsize=14)
    ax[1].set_xlabel("time (min)", fontsize=14)
    ax[0].legend()
    for x in ax:
        for tick in x.xaxis.get_major_ticks():
             tick.label.set_fontsize(14) 
        for tick in x.yaxis.get_major_ticks():
             tick.label.set_fontsize(14) 

    fig.savefig("Ca_and_Fluo_for_Glu_NMDA_bath_no_ER.png", dpi=100)
    plt.show()
        

    

