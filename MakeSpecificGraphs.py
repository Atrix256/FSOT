import sys
from PIL import Image
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from scipy import signal

fileNames = glob.glob("out/projective_*.csv")

for fileName in fileNames:

    if ".movement." in fileName:
        continue
    
    print(fileName)

    # get the values
    df = pd.read_csv(fileName)
    valuesX = df["X"].to_numpy()
    valuesY = df["Y"].to_numpy()

    # plot them
    pointsX = np.zeros(valuesX.shape[0] * 10)
    pointsY = np.zeros(valuesY.shape[0] * 10)

    for v in valuesX:
        pointsX[min(max(int((v * 0.5 + 0.5) * pointsX.shape[0]), 0), pointsX.shape[0]-1)] = 1.0

    for v in valuesY:
        pointsY[min(max(int((v * 0.5 + 0.5) * pointsY.shape[0]), 0), pointsY.shape[0]-1)] = 1.0

    # DFT the plots and zero out DC
    dftX = abs(np.fft.rfft(pointsX))
    dftX[0] = 0
    #dftX = np.fft.fftshift(dftX)

    dftY = abs(np.fft.rfft(pointsY))
    dftY[0] = 0
    #dftY = np.fft.fftshift(dftY)

    # normalize the amplitudes
    #maxAmplitude = max(np.max(dftX), np.max(dftY))
    #dftX /= maxAmplitude
    #dftY /= maxAmplitude

    # Graph the DFTs

    dftX = signal.resample(dftX, 100)
    dftY = signal.resample(dftY, 100)
    t = np.arange(0, 1, 1.0 / dftX.shape[0])
    
    fig = plt.figure()
    plt.xlim(0.0, 1.0)
    #plt.ylim(0.0, 1.0)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Amplitude')
    plt.plot(t, dftX, label = "X")
    plt.plot(t, dftY, label = "Y")
    #plt.yticks([])

    plt.legend()

    outFileName =  os.path.splitext(fileName)[0] + ".magnitudeXY.png"
    
    #plt.show()

    fig.tight_layout()
    fig.savefig(outFileName, bbox_inches='tight')

'''

fileNames =[
    "out/batch1_square.csv",
    "out/batch4_square.csv",
    "out/batch16_square.csv",
    "out/batch64_square.csv",
    "out/batch256_square.csv",
]

fig, ax = plt.subplots()
plt.title('Log/Log Average Movement Each Iteration')

xValues = range(1, 64000, 634)

for fileName in fileNames:
    print(fileName)
    df = pd.read_csv(fileName).drop(['Iteration'], axis=1)
    ax.plot(xValues, df['Avg. Movement'], label=fileName)

ax.legend()

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig("out/batch_square.graph.png", bbox_inches='tight')



fileNames =[
    "out/batch1_stratified_square.csv",
    "out/batch4_stratified_square.csv",
    "out/batch16_stratified_square.csv",
    "out/batch64_stratified_square.csv",
    "out/batch256_stratified_square.csv",
]

fig, ax = plt.subplots()
plt.title('Log/Log Average Movement Each Iteration')

xValues = range(1, 64000, 634)

for fileName in fileNames:
    print(fileName)
    df = pd.read_csv(fileName).drop(['Iteration'], axis=1)
    ax.plot(xValues, df['Avg. Movement'], label=fileName)

ax.legend()

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig("out/batch_stratified_square.graph.png", bbox_inches='tight')



fileNames =[
    "out/square.csv",
    "out/circle.csv",
    "out/square_GR.csv",
    "out/circle_GR.csv",
]

fig, ax = plt.subplots()
plt.title('Log/Log Average Movement Each Iteration')

xValues = range(1, 64000, 634)

for fileName in fileNames:
    print(fileName)
    df = pd.read_csv(fileName).drop(['Iteration'], axis=1)
    ax.plot(xValues, df['Avg. Movement'], label=fileName)

ax.legend()

fig.axes[0].set_xscale('log', base=2)
fig.axes[0].set_yscale('log', base=2)

fig.tight_layout()
fig.savefig("out/GR.graph.png", bbox_inches='tight')
'''
