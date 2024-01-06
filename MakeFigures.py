import sys
from PIL import Image
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

makeTiles = False

fileNames = glob.glob("out/*.points.png") + glob.glob("out/*.gauss.png") + glob.glob("out/*.texture2D.png")

# Make DFTs of all images
for fileName in fileNames:

    print(fileName)

    fileNameNoExtension = os.path.splitext(fileName)[0]

    # Load the image
    loadedImage = Image.open(fileName)

    # DFT
    if ".gauss." not in fileName:
        # convert to greyscale to handle color images
        im = np.array(loadedImage.convert('L'), dtype=float) / 255.0

        # convert point images to black and white. we could convert to mode '1' but that does dithering and we don't want that
        if ".texture2D." not in fileName:
            im2 = (im == 1.0) * 1.0
            im = im2

        # get the DFT magnitude, zero out DC and shift it so DC is in the middle
        dft = abs(np.fft.fft2(im))
        dft[0,0] = 0.0
        dft = np.fft.fftshift(dft)

        # log and normalize DFT
        imOut = np.log(1+dft)
        themin = np.min(imOut)
        themax = np.max(imOut)
        if themin != themax:
            imOut = (imOut - themin) / (themax - themin)
        else:
            imOut = imOut

        # Write out DFT
        outFileName = fileNameNoExtension + ".magnitude.png"
        Image.fromarray((imOut*255.0).astype(np.uint8), mode="L").save(outFileName)
    
    # Tile image as 3x3
    if makeTiles:
        imOut = Image.new(loadedImage.mode, (loadedImage.size[0] * 3, loadedImage.size[1] * 3), 255)
        for i in range(3):
            for j in range(3):
                imOut.paste(loadedImage, (i*loadedImage.size[0],j*loadedImage.size[1]))
        imOut.save(fileNameNoExtension + ".3x3.png")
    
    # Tile image as 11x11
    if makeTiles:
        imOut = Image.new(loadedImage.mode, (loadedImage.size[0] * 11, loadedImage.size[1] * 11), 255)
        for i in range(11):
            for j in range(11):
                imOut.paste(loadedImage, (i*loadedImage.size[0],j*loadedImage.size[1]))
        imOut.save(fileNameNoExtension + ".11x11.png")


fileNames = glob.glob("out/*.movement.csv")

for fileName in fileNames:
    print(fileName)

    fileNameNoExtension = os.path.splitext(fileName)[0]

    fig, ax = plt.subplots()
    df = pd.read_csv(fileName).drop(['Iteration'], axis=1)

    ax.plot(df['Avg. Movement'], label="Avg. Movement")

    plt.title('Log/Log Average Movement Each Iteration: ' + fileName)

    fig.axes[0].set_xscale('log', base=2)
    fig.axes[0].set_yscale('log', base=2)

    fig.tight_layout()
    fig.savefig(fileNameNoExtension + ".graph.png", bbox_inches='tight')
