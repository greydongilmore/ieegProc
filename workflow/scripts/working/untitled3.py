#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 17:13:36 2025

@author: greydon
"""
import binascii
import numpy as np
import sys,zlib


def hasher(fileName):
	f = open(ros_fname,"rb")
	crc = 0
	while True:
		buffer = f.read(1024*1024)
		if len(buffer) == 0:
			f.close()
			return crc
		crc = binascii.crc32(buffer, crc)

ros_fname=r"/home/greydon/Documents/data/emory_seeg/derivatives/slicer_scene/sub-EMOP0340/Maynard_Kristofer_714975_2022-08-30_ROSA/MAYNARD KRISTOFER    20220830 080609.ros"




def crc32(v):
	return binascii.crc32(v.encode("utf8"))

crc = hasher(ros_fname)

with open(ros_fname, 'rb') as f:
	textfile = f.read()

print(binascii.crc32(textfile))


asciiText = removeNonAscii(textfile)
for ib in textfile:
	crc = binascii.crc32(ib, crc)

iHash=hasher(ros_fname)
sHash = '%08X' % iHash
crc_max = int("FFFFFFFF", base=16)
crc_min = int("00000000", base=16)
targetCRC = int(str(crc), base=16)

newContents = calcNewContents(targetCRC, crc)
print('Four bytes to append (hex): %s' % itos(newContents))
