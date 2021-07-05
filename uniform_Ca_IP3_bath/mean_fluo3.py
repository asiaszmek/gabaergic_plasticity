import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+',
                        help='input files')
parser.add_argument("--stim_start", help="Stimulus onset (ms)",
                    type=float)
specie = "Fluo3Ca"

if __name__ == "__main__":

    data = []
    args = parser.parse_args()
    files = []
    for name in args.input:
        files.append(name)
    if not files:
        sys.exit('Do specify at least one totals filename')
    t_stim = args.stim_start
    for fname in files:
        f = open(fname, "r")
        header = f.readline().split()
        datas = np.loadtxt(f)
        time = datas[:, 0]
        idx = header.index(specie)
        data.append(datas[:, idx])
    max_len = max([len(d) for d in data])
    out = np.zeros((max_len, len(data)))
    for i, d in enumerate(data):
        out[:len(d), i] = np.array(d)
    mean_signal = out.mean(axis=1)
    dt = time[1] - time[0]

    mean_fluo = mean_signal[0:int(t_stim//dt)].mean()

    fig, ax = plt.subplots(1, 1)
    ax.plot(time/1000, (mean_signal-mean_fluo)/mean_fluo)
    ax.set_xlabel("time (sec)")
    ax.set_ylabel("% Fluorescence")
    plt.show()
    
        
                   

