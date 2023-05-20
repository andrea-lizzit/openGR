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
import json
from scipy import ndimage
import scipy
from tqdm import tqdm

def sph2pix(map, width, height):
	return map * np.array([height/np.pi, width / (2 * np.pi), 1])

with open(args.map, 'r') as f:
	map = np.array(json.load(f))

# load images
ucs = Image.open(args.ucs)
lcs = Image.open(args.lcs)
# convert PIL image to numpy array
aucs = np.array(ucs)
alcs = np.array(lcs)

def upsample(map, width, height, closed=False):
	if closed:
		width += 1
		height += 1
	# upsample pixmap using map_coordinates
	finegrid = np.array(np.meshgrid(np.arange(height), np.arange(width), np.arange(3), indexing='ij'))
	upsample_map = finegrid * np.reshape(np.array([(map.shape[0] - 1) / (args.height-1), (map.shape[1] - 1) / (args.width-1), 1]), (3, 1, 1, 1))
	upsampled_pixmap = ndimage.map_coordinates(map, upsample_map, order=1)
	if closed:
		upsampled_pixmap = upsampled_pixmap[:, :-1, :-1]
	return upsampled_pixmap


def warp_image(ucs, lcs, map):
	# convert spherical coordinates to pixel coordinates
	ucs_pixmap = sph2pix(map, ucs.shape[1], ucs.shape[0])
	# if args.width is None:
	# 	args.width = ucs.shape[1]
	# if args.height is None:
	# 	args.height = ucs.shape[0]
	if args.height is not None and args.width is not None:
		ucs_pixmap = upsample(ucs_pixmap, args.width, args.height)
	ucs_pixmap = ucs_pixmap.transpose((2, 0, 1))
	lcs_pixmap = sph2pix(map, lcs.shape[1], lcs.shape[0])
	# if args.width is None:
	# 	args.width = lcs.shape[1]
	# if args.height is None:
	# 	args.height = lcs.shape[0]
	if args.height is not None and args.width is not None:
		lcs_pixmap = upsample(lcs_pixmap, args.width, args.height)
	lcs_pixmap = lcs_pixmap.transpose((2, 0, 1))
	lcs_pixmap = np.flip(lcs_pixmap, axis=2)
	ucs_channels, lcs_channels = np.split(ucs, 3, axis=2), np.split(lcs, 3, axis=2)
	warped_channels = []
	for uc, lc in zip(ucs_channels, lcs_channels):
		uc, lc = np.squeeze(uc), np.squeeze(lc)
		warped_uc = ndimage.map_coordinates(uc, ucs_pixmap[:2], order=1)
		warped_lc = ndimage.map_coordinates(lc, lcs_pixmap[:2], order=1)
		mask = np.rint(ucs_pixmap[2]).astype(np.uint8)
		warped_skies = warped_uc * mask + warped_lc * (1 - mask)
		warped_channels.append(warped_skies)
	warped_image = np.stack(warped_channels, axis=2)
	return warped_image

warped_image = warp_image(aucs, alcs, map)
Image.fromarray(warped_image, mode='RGB').save(args.output)