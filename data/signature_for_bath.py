from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

colors = {
    "100 uM Glu": "tab:orange",
    "50 uM NMDA": "tab:blue",
    "10 uM Glu": "tab:green",
}

fnames_spine = {
    "10 uM Glu": [
        "model_150_nM_spine_Glu_10_uM_trial0_sa1[0].txt",
        "model_150_nM_spine_Glu_10_uM_trial1_sa1[0].txt",
        "model_150_nM_spine_Glu_10_uM_trial2_sa1[0].txt",
        "model_150_nM_spine_Glu_10_uM_trial3_sa1[0].txt",
    ],
    "50 uM NMDA": [
        "model_150_nM_spine_NMDA_trial0_sa1[0].txt",
        "model_150_nM_spine_NMDA_trial1_sa1[0].txt",
        "model_150_nM_spine_NMDA_trial2_sa1[0].txt",
        "model_150_nM_spine_NMDA_trial3_sa1[0].txt",
    ],
    "100 uM Glu": [
        "model_250_nM_spine_Glu_100_uM_trial0_sa1[0].txt",
        "model_250_nM_spine_Glu_100_uM_trial1_sa1[0].txt",
        "model_250_nM_spine_Glu_100_uM_trial2_sa1[0].txt",
        "model_250_nM_spine_Glu_100_uM_trial3_sa1[0].txt",

    ],
}


fnames_dend = {
    "10 uM Glu": [
        "model_150_nM_spine_Glu_10_uM_trial0_dend_.txt",
        "model_150_nM_spine_Glu_10_uM_trial1_dend_.txt",
        "model_150_nM_spine_Glu_10_uM_trial2_dend_.txt",
        "model_150_nM_spine_Glu_10_uM_trial3_dend_.txt",
    ],
    "50 uM NMDA": [
        "model_150_nM_spine_NMDA_trial0_dend_.txt",
        "model_150_nM_spine_NMDA_trial1_dend_.txt",
        "model_150_nM_spine_NMDA_trial2_dend_.txt",
        "model_150_nM_spine_NMDA_trial3_dend_.txt",
    ],
    "100 uM Glu": [
        "model_250_nM_spine_Glu_100_uM_trial0_dend_.txt",
        "model_250_nM_spine_Glu_100_uM_trial1_dend_.txt",
        "model_250_nM_spine_Glu_100_uM_trial2_dend_.txt",
        "model_250_nM_spine_Glu_100_uM_trial3_dend_.txt",

    ],
}






data_dend = OrderedDict()
data_spine = OrderedDict()
keys = ["10 uM Glu", "50 uM NMDA", "100 uM Glu"]

def make_figs(dend, spine, idexes, ylabel, fname):
    new_dend = {}
    new_spine = {}
    time_idx = 0
    fig, ax = plt.subplots(2, len(keys), figsize= (24, 12))
    i = 0 # dendrite
    min_val = 100000
    max_val = 0
    for i, key in enumerate(keys):
        avg = []
        new_dend[key] = []
        for j, d in enumerate(dend[key]):
            time = (dend[key][j][:, 0]-100000)/60000 #convert to min stim
            out = np.zeros((len(time)))
            for idx in idexes:
                out += dend[key][j][:, idx]
          
            if out.max() > max_val:
                max_val = out.max()
            if out.min() < min_val:
                min_val = out.min()
            ax[0][i].plot(time, out, label="trial %d" % j)
            avg.append(out)
            new_dend[key].append(out)
        out_avg = np.array(avg).mean(axis=0)
        #ax[0][i].plot(time, out_avg, label="mean CaMKII")
        ax[0][i].set_title("%s dendrite" % key, fontsize=14)
        if i:
            ax[0][i].set_yticks([])
    for x in ax[0]:
        x.set_ylim([min_val , max_val])
        x.set_xticks([])
        
    for i, key in enumerate(keys):
        max_val = 0
        min_val = 1000000
        new_spine[key] = []
        for j, d in enumerate(spine[key]):
            time = (dend[key][j][:, 0]-100000)/60000 #convert to min stim
            out = np.zeros((len(time)))
            for idx in idexes:
                out += spine[key][j][:, idx]
            if out.max() > max_val:
                max_val = out.max()
            if out.min() < min_val:
                min_val = out.min()
            ax[1][i].plot(time, out, label="trial %d" % j)
            avg.append(out)
            new_spine[key].append(out)
        out_avg = np.array(avg).mean(axis=0)
        #ax[1][i].plot(time, out_avg, label="mean CaMKII")
        ax[1][i].set_title("%s spine" % key, fontsize=14)
        if i:
            ax[1][i].set_yticks([])

    for x in ax[1]:
        x.set_ylim([min_val * 0.95, max_val*1.05])
    ax[0][0].set_ylabel(ylabel, fontsize=14)
    ax[1][0].set_ylabel(ylabel, fontsize=14)
    ax[1][0].set_xlabel("time (min)", fontsize=14)
    for x in ax[0]:
        for tick in x.xaxis.get_major_ticks():
             tick.label.set_fontsize(14) 
        for tick in x.yaxis.get_major_ticks():
             tick.label.set_fontsize(14)
    for x in ax[1]:
        for tick in x.xaxis.get_major_ticks():
             tick.label.set_fontsize(14) 
        for tick in x.yaxis.get_major_ticks():
             tick.label.set_fontsize(14) 
 
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    return new_dend, new_spine

if __name__ == "__main__":
    header_old = []
    for key in keys:
        data_spine[key] = []
        data_dend[key] = []
        for fname in fnames_spine[key]:
            f = open(fname)
            header = f.readline().split()
            data = np.loadtxt(f)
            time = data[:, 0]
            data_spine[key].append(data)
            if header != header_old:
                header_old = header
                print("Header has changed")
        for fname in fnames_dend[key]:
            f = open(fname)
            header = f.readline().split()
            data = np.loadtxt(f)
            data_dend[key].append(data)
            if header != header_old:
                header_old = header
                print("Header has changed")

    all_CK_act = []
    for i, specie in enumerate(header):
        if "CKC" in specie or "CKp" in specie:
            all_CK_act.append(i)
    print(all_CK_act)
    CK_act = []
    for i, specie in enumerate(header):
        if "CKC" in specie:
            CK_act.append(i)
    CKp = []
    for i, specie in enumerate(header):
        if "CKp" in specie:
            CKp.append(i)

    PP2B_all = []
    for i, specie in enumerate(header):
        if "PP2BCaMC" in specie:
            PP2B_all.append(i)
    
        if "PP2B" in specie and "GluR" in specie:
            PP2B_all.append(i)

    PP2B = []
    for i, specie in enumerate(header):
        if "PP2BCaMCa4" in specie:
            PP2B.append(i)
    
        if "PP2B" in specie and "GluR" in specie:
            PP2B.append(i)
    
    PP1 = []
    for i, specie in enumerate(header):
        if specie == "PP1":
            PP1.append(i)
        if "PP1" in specie and "GluR" in specie:
            PP1.append(i)
        if "PP1" in specie and "CK" in specie:
            PP1.append(i)

    epac = []
    for i, specie in enumerate(header):
        if specie == "Epac1cAMP":
            epac.append(i)

    acamkii_dend, acamkii_spine = make_figs(data_dend, data_spine,
                                            all_CK_act, "Active CaMKII (nM)",
                                            "Bath_paradigms_active_CaMKII.png")
    pcamkii_dend, pcamkii_spine = make_figs(data_dend, data_spine,
                                            CKp, "pCaMKII (nM)",
                                            "Bath_paradigms_pCaMKII.png")
    camkii_dend, camkii_spine = make_figs(data_dend, data_spine,
                                          CK_act, "CaMCa4 bound CaMKII (nM)",
                                            "Bath_paradigms_CaMKII.png")
    pp2ball_dend, pp2ball_spine = make_figs(data_dend, data_spine,
                                            PP2B_all, "PP2B bound CaM (nM)",
                                            "Bath_paradigms_PP2B_bound.png")
    pp2b_dend, pp2b_spine = make_figs(data_dend, data_spine,
                                      PP2B, "Activated PP2B (nM)",
                                            "Bath_paradigms_PP2B.png")
    pp1_dend, pp1_spine = make_figs(data_dend, data_spine,
                                    PP1, "PP1 (nM)",
                                    "Bath_paradigms_PP1.png")
    
    epac_dend, epac_spine = make_figs(data_dend, data_spine,
                                      epac, "Epac1cAMP (nM)",
                                      "Bath_paradigms_epac.png")
    
    
    fig, ax = plt.subplots(1, 2)
    min_val = 100000
    max_val = 0
    for key in keys:
        avg_pcamkii_dend = np.array(pcamkii_dend[key]).mean(axis=0)
        avg_pcamkii_spine = np.array(pcamkii_spine[key]).mean(axis=0)
        avg_acamkii_dend = np.array(acamkii_dend[key]).mean(axis=0)
        avg_acamkii_spine = np.array(acamkii_spine[key]).mean(axis=0)
        avg_camkii_dend = np.array(camkii_dend[key]).mean(axis=0)
        avg_camkii_spine = np.array(camkii_spine[key]).mean(axis=0)
        avg_pp2ball_dend = np.array(pp2ball_dend[key]).mean(axis=0)
        avg_pp2ball_spine = np.array(pp2ball_spine[key]).mean(axis=0)
        avg_pp2b_dend = np.array(pp2b_dend[key]).mean(axis=0)
        avg_pp2b_spine = np.array(pp2b_spine[key]).mean(axis=0)
        avg_pp1_dend = np.array(pp1_dend[key]).mean(axis=0)
        avg_pp1_spine = np.array(pp1_spine[key]).mean(axis=0)
        avg_epac_dend = np.array(epac_dend[key]).mean(axis=0)
        avg_epac_spine = np.array(epac_spine[key]).mean(axis=0)
        dend = avg_acamkii_dend/(avg_pp2ball_dend+avg_epac_dend)
        spine = avg_acamkii_spine/(avg_pp2ball_spine+avg_epac_spine)
        ax[1].plot((time-100000)/60000, dend,
                   label=key, color=colors[key])
        ax[0].plot((time-100000)/60000, spine,
                   label=key, color=colors[key])
        if dend.max() > max_val:
            max_val = dend.max()
        if spine.max() > max_val:
            max_val = spine.max()
        if dend.min() < min_val:
            min_val = dend.min()
        if spine.min() < min_val:
            min_val = spine.min()
    ax[0].legend()
    ax[0].set_xlabel("time (min)", fontsize=14)
    ax[0].set_ylabel("active CaMKII/(active PP2B + Epac)", fontsize=14)
    ax[1].set_title("dendrite", fontsize=14)
    ax[1].set_yticklabels([])
    ax[0].set_title("spine", fontsize=14)
    for x in ax:
        for tick in x.xaxis.get_major_ticks():
             tick.label.set_fontsize(14) 
        for tick in x.yaxis.get_major_ticks():
             tick.label.set_fontsize(14) 
        x.set_ylim([min_val *0.95, 1.05*max_val])
    fig.savefig("CaMKII_to_PP2B_and_Epac.png", bbox_inches="tight")
    plt.show()

    
