from enum import Enum

import requests
import straw
import colorsys
from PIL import Image
import numpy as np
import pandas as pd
import math
from pandarallel import pandarallel
import functools
from copy import deepcopy as dc
import scipy.ndimage
from skimage import feature


def CreateMapFromHic(
		HicFile: str,
		ChromX: str,
		ChromY: str,
		Normalization: str,
		BinSize: int,
		StartX: int = None,
		EndX: int = None,
		StartY: int = None,
		EndY: int = None,
		) -> pd.DataFrame:
	
	NormTypes = ["NONE", "VC", "VC_SQRT", "KR"]
	BinTypes = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 500, 200, 100, 50, 20, 5, 2, 1]
	if Normalization not in NormTypes: raise ValueError(f"Wrong normalization mode: '{str(Normalization)}'. Valid values: {', '.join([str(i) for i in NormTypes])}")
	if BinSize not in BinTypes: raise ValueError(f"Wrong bin size: '{str(BinSize)}'. Valid values: {', '.join([str(i) for i in BinTypes])}")
	if (StartX is not None) and (EndX is not None) and ((StartX >= EndX) or (StartX <= 0) or (EndX <= 0)): raise ValueError(f"Wrong coordinates: {StartX=}, {EndX=}")
	if (StartY is not None) and (EndY is not None) and ((StartY >= EndY) or (StartY <= 0) or (EndY <= 0)): raise ValueError(f"Wrong coordinates: {StartY=}, {EndY=}")
	
	BinMode = "BP" if BinSize > 1000 else "FRAG"
	CoordsX = ChromX + ("" if ((StartX is None) or (EndX is None)) else f":{str(StartX)}:{str(EndX)}")
	CoordsY = ChromY + ("" if ((StartY is None) or (EndY is None)) else f":{str(StartY)}:{str(EndY)}")
	Straw = straw.straw(Normalization, HicFile, CoordsX, CoordsY, BinMode, BinSize)
	Data = pd.DataFrame(Straw).transpose().set_index([0, 1]).unstack().fillna(0.0).applymap(int).rename_axis(None, axis=0).sort_index()
	Data.columns = sorted([int(item[1]) for item in Data.columns])
	
	return Data

def MaxContacts(Data: pd.DataFrame) -> int:
	Max = Data.max().max()
	return Max if Max != 0 else 1

def SaveMapToPngRainbow(
		Dataframe: pd.DataFrame,
		PngFileName: str) -> None:
	Max = MaxContacts(Dataframe)
	Dataframe = Dataframe.applymap(lambda x: list(colorsys.hsv_to_rgb(1.0 - (x / Max), 0.8, x / Max)))
	Dataframe = np.array(Dataframe.applymap(lambda x: [np.byte(int(item * 255)) for item in x]).values.tolist())
	Image.fromarray(Dataframe, 'RGB').save(PngFileName)
	
def SaveMapToPng3Channels(
		Dataframes: dict,
		PngFileName: str) -> None:
	Max = {index: MaxContacts(item) for index, item in Dataframes.items()}
	Dataframes = {index: item.applymap(lambda x: [np.byte(int(x / Max[index] * 255))]) for index, item in Dataframes.items()}
	Data = np.array(Dataframes['R'].add(Dataframes['G']).add(Dataframes['B']).values.tolist())
	Image.fromarray(Data, 'RGB').save(PngFileName)




MapR = CreateMapFromHic(
		HicFile = "/home/fairwind/GSE168470_KARPAS_A485_merged_rep1234.hic",
		ChromX= "1",
		ChromY="1",
		Normalization= "VC_SQRT",
		BinSize= 500000
		)

MapG = MapR.applymap(lambda x: 0)
MapB = MapR.applymap(lambda x: 0)

MapR = MapR.applymap(lambda x: math.log(x + 1))

SaveMapToPng3Channels({'R': MapR, 'G': MapG, 'B': MapB}, 'test.png')

