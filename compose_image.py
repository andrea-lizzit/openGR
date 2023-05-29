import argparse

parser = argparse.ArgumentParser()
parser.add_argument("map", type=str, help='Coordinate map')
parser.add_argument("ucs", type=str, help='Image of the upper celestial sphere')
parser.add_argument("lcs", type=str, help='Image of the lower celestial sphere')
parser.add_argument("output", type=str, help='Output image')
parser.add_argument("--width", type=int, default=None, help='Width of the output image')
parser.add_argument("--height", type=int, default=None, help='Height of the output image')

args = parser.parse_args()

from PIL import Image
import numpy as np
import h5py
from scipy import ndimage
import scipy
from tqdm import tqdm

def sph2pix(map, width, height):
	closed_pixmap = map * np.array([(height)/np.pi, (width) / (2 * np.pi), 1])
	return closed_pixmap[:-1, :-1, :]

def differential(map):
	dx = np.diff(np.concatenate((map, map[:, :1]), axis=1), axis=1)
	dy = np.diff(np.concatenate((map, map[:1, :]), axis=0), axis=0)
	return dy, dx

f = h5py.File(args.map, 'r')
map = np.array(f['map'])

# load images
ucs = Image.open(args.ucs)
lcs = Image.open(args.lcs)
# convert PIL image to numpy array
aucs = np.array(ucs)
alcs = np.array(lcs)

def upsample(map, width, height):
	# upsample pixmap using map_coordinates
	finegrid = np.array(np.meshgrid(np.arange(height), np.arange(width), np.arange(3), indexing='ij'))
	upsample_map = finegrid * np.reshape(np.array([(map.shape[0] - 1) / (height-1), (map.shape[1] - 1) / (width-1), 1]), (3, 1, 1, 1))
	upsampled_pixmap = ndimage.map_coordinates(map, upsample_map, order=1)
	return upsampled_pixmap


def warp_image(ucs, lcs, map):
	# convert spherical coordinates to pixel coordinates
	if args.height is not None and args.width is not None:
		map = upsample(map, args.width+1, args.height+1) # +1 because the edges will be discarded
	ucs_pixmap = sph2pix(map, ucs.shape[1], ucs.shape[0])
	lcs_pixmap = sph2pix(map, lcs.shape[1], lcs.shape[0])
	ucs_pixmap = (ucs_pixmap.transpose((2, 0, 1)))
	lcs_pixmap = (lcs_pixmap.transpose((2, 0, 1)))
	lcs_pixmap = np.flip(lcs_pixmap, axis=2)
	ucs_channels, lcs_channels = np.split(ucs, 3, axis=2), np.split(lcs, 3, axis=2)
	warped_channels = []
	for uc, lc in zip(ucs_channels, lcs_channels):
		uc, lc = np.squeeze(uc), np.squeeze(lc)
		warped_uc = ndimage.map_coordinates(uc, ucs_pixmap[:2], order=1, cval=127, mode="grid-wrap")
		warped_lc = ndimage.map_coordinates(lc, lcs_pixmap[:2], order=1, cval=127, mode="grid-wrap")
		mask = np.rint(ucs_pixmap[2]).astype(np.uint8)
		warped_skies = warped_uc * mask + warped_lc * (1 - mask)
		warped_channels.append(warped_skies)
	warped_image = np.stack(warped_channels, axis=2)
	return warped_image

warped_image = warp_image(aucs, alcs, map)
Image.fromarray(warped_image, mode='RGB').save(args.output)