import collections
import pprint

import pytest

import numpy as np
import qcelemental as qcel
from qcelemental.testing import compare, compare_values

from .utils import using

import qcdb


data = collections.defaultdict(dict)

data["C2H2"]["ref"] = [
    [0.000000000000, 0.000000000000, 0.650000000000],
    [0.000000000000, 0.000000000000, -0.650000000000],
    [0.000000000000, 0.000000000000, 1.750000000000],
    [0.000000000000, 0.000000000000, -1.750000000000],
]


data["N2"]["ref"] = [
    [0.000000000000, 0.000000000000, -0.550000000000],
    [0.000000000000, 0.000000000000, 0.550000000000],
]


data["CN"]["ref"] = [
    [0.000000000000, 0.000000000000, -0.753922540200],
    [0.000000000000, 0.000000000000, 0.646077459800],
]


data["HCCCl"]["ref"] = [
    [0.000000000000, 0.000000000000, -2.719188122301],
    [0.000000000000, 0.000000000000, -1.719188122301],
    [0.000000000000, 0.000000000000, -0.619188122301],
    [0.000000000000, 0.000000000000, 0.880811877699],
]

data["CHFClBr"]["ref"] = [
    [0.251995331546, -0.251226142375, -1.228834841943],
    [0.251995331546, -0.251226142375, -0.228834841943],
    [0.251995331546, -1.217151968664, 0.029984203160],
    [1.088511635284, 0.231736770769, 0.029984203160],
    [-0.584520972191, 0.231736770769, 0.029984203160],
]


data["CH2ClBr"]["ref"] = [
    [-0.588382367674, 0.890373072017, 0.000000000000],
    [-0.588382367674, -0.109626927983, 0.000000000000],
    [0.377543458615, -0.368445973085, 0.000000000000],
    [-1.071345280819, -0.368445973085, 0.836516303738],
    [-1.071345280819, -0.368445973085, -0.836516303738],
]


data["HOCl"]["ref"] = [
    [-1.074855537230, 1.371823577282, 0.000000000000],
    [-1.074855537230, 0.371823577282, 0.000000000000],
    [0.522621918106, -0.209610666372, 0.000000000000],
]


data["C4H4Cl2F2"]["ref"] = [
    [0.432781050498, 1.898774028282, 0.810337938486],
    [-1.658744642774, 0.805191018766, -0.984829058337],
    [1.658744642774, -0.805191018766, 0.984829058337],
    [-0.432781050498, -1.898774028282, -0.810337938486],
    [-0.317971784026, 2.532165941971, 2.640915161238],
    [-1.615729990528, 1.614062700629, -2.881498569657],
    [1.615729990528, -1.614062700629, 2.881498569657],
    [0.317971784026, -2.532165941971, -2.640915161238],
    [-4.852178875691, 1.024620478757, 0.190249941464],
    [4.852178875691, -1.024620478757, -0.190249941464],
    [-1.913713787211, -3.739567959534, 0.258534542158],
    [1.913713787211, 3.739567959534, -0.258534542158],
]


data["HOOH_dimer"]["ref"] = [
    [0.991126228500, -1.797922633300, 0.146518251500],
    [2.769109309500, -1.348521864900, -0.007155768400],
    [2.517803031100, 1.380837492300, -0.115405801400],
    [3.288320045300, 1.830859509500, 1.475770682500],
    [-3.288320045300, -1.830859509500, -1.475770682500],
    [-2.517803031100, -1.380837492300, 0.115405801400],
    [-2.769109309500, 1.348521864900, 0.007155768400],
    [-0.991126228500, 1.797922633300, -0.146518251500],
]


data["HOOH"]["ref"] = [
    [-0.657774170895, 0.990250821217, -0.765560415408],
    [-0.134283313545, 0.737880743551, 0.048237265941],
    [0.134283313545, -0.737880743551, 0.048237265941],
    [0.657774170895, -0.990250821217, -0.765560415408],
]


data["NOHOHOH"]["ref"] = [
    [0.000000000000, 0.000000000000, -1.219526050935],  # Gh
    [0.000000000000, 0.000000000000, -0.219526050935],
    [0.000000000000, -1.477211629518, 0.040946215565],
    [1.279302797929, 0.738605814759, 0.040946215565],
    [-1.279302797929, 0.738605814759, 0.040946215565],
    [0.899871989900, -1.767037453836, 0.366877793280],
    [1.080363329511, 1.662830730326, 0.366877793280],
    [-1.980235319411, 0.104206723511, 0.366877793280],
]


data["H2O"]["ref"] = [
    [0.000000000000, 0.816641555162, -0.512554059234],
    [0.000000000000, 0.000000000000, 0.064591130803],
    [0.000000000000, -0.816641555162, -0.512554059234],
]


data["CH2F2"]["ref"] = [
    [0.000000000000, 0.000000000000, 1.089095845660],
    [0.000000000000, -2.122315581200, -0.459816147540],
    [-0.000000000000, 2.122315581200, -0.459816147540],
    [1.708413985000, 0.000000000000, 2.184106800160],
    [-1.708413985000, -0.000000000000, 2.184106800160],
]


data["NH3"]["ref"] = [
    [0.000000000000, -1.071293777318, 0.000000000000],  # Gh
    [0.000000000000, -0.071293777318, 0.000000000000],
    [-0.430496198842, 0.330193571336, -0.745641288860],
    [-0.430496198842, 0.330193571336, 0.745641288860],
    [0.860992397685, 0.330193571336, 0.000000000000],
]


data["BrF5"]["ref"] = [
    [0.000000000000, 0.000000000000, -1.514287735691],
    [0.000000000000, 0.000000000000, 0.185712264309],
    [0.000000000000, -1.700000000000, 0.185712264309],
    [-1.700000000000, 0.000000000000, 0.185712264309],
    [1.700000000000, 0.000000000000, 0.185712264309],
    [0.000000000000, 1.700000000000, 0.185712264309],
]


data["N2H2"]["ref"] = [
    [0.000000000000, 0.700000000000, 0.000000000000],
    [0.000000000000, -0.700000000000, 0.000000000000],
    [-0.642787609687, 1.466044443119, 0.000000000000],
    [0.642787609687, -1.466044443119, 0.000000000000],
]


data["NOHOHOH"]["ref"] = [
    [0.000000000000, 0.000000000000, -1.000000000000],  # Gh
    [0.000000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, -1.500000000000, 0.000000000000],
    [1.299038105677, 0.750000000000, 0.000000000000],
    [-1.299038105677, 0.750000000000, 0.000000000000],
    [0.939692620786, -1.842020143326, 0.000000000000],
    [1.125389928010, 1.734807753012, 0.000000000000],
    [-2.065082548796, 0.107212390313, 0.000000000000],
]


data["TFCOT"]["ref"] = [
    [-1.618188000000, -0.437140000000, -0.409373000000],
    [-1.394411000000, 0.896360000000, -0.429596000000],
    [-0.896360000000, -1.394411000000, 0.429596000000],
    [-0.437140000000, 1.618188000000, 0.409373000000],
    [0.437140000000, -1.618188000000, 0.409373000000],
    [0.896360000000, 1.394411000000, 0.429596000000],
    [1.394411000000, -0.896360000000, -0.429596000000],
    [1.618188000000, 0.437140000000, -0.409373000000],
    [2.147277000000, -1.690111000000, -1.235043000000],
    [1.690111000000, 2.147277000000, 1.235043000000],
    [-2.147277000000, 1.690111000000, -1.235043000000],
    [-1.690111000000, -2.147277000000, 1.235043000000],
    [0.878010000000, -2.418132000000, 1.029595000000],
    [-2.418132000000, -0.878010000000, -1.029595000000],
    [-0.878010000000, 2.418132000000, 1.029595000000],
    [2.418132000000, 0.878010000000, -1.029595000000],
]


data["Li_H2O_4_p"]["ref"] = [
    [0.000000000000, 0.000000000000, -1.000000000000],  # Gh
    [0.000000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, -1.000000000000, 0.000000000000],  # Gh
    [0.000000000000, 0.000000000000, 1.000000000000],  # Gh
    [-1.497220431853, 0.000000000000, -1.169756803119],
    [1.497220431853, 0.000000000000, -1.169756803119],
    [0.000000000000, -1.497220431853, 1.169756803119],
    [0.000000000000, 1.497220431853, 1.169756803119],
    [-1.565808146965, 0.498804253130, -1.977713511362],
    [-2.264066335924, -0.551608165437, -1.068881676192],
    [1.565808146965, -0.498804253130, -1.977713511362],
    [2.264066335924, 0.551608165437, -1.068881676192],
    [-0.498804253130, -1.565808146965, 1.977713511362],
    [0.551608165437, -2.264066335924, 1.068881676192],
    [0.498804253130, 1.565808146965, 1.977713511362],
    [-0.551608165437, 2.264066335924, 1.068881676192],
]


data["ethylene_cation"]["ref"] = [
    [0.000000000000, 0.000000000000, -0.705000000000],
    [0.000000000000, 0.000000000000, 0.705000000000],
    [0.353742012310, 0.854008763700, -1.282611998014],
    [-0.353742012310, -0.854008763700, -1.282611998014],
    [-0.353742012310, 0.854008763700, 1.282611998014],
    [0.353742012310, -0.854008763700, 1.282611998014],
]


data["ethane_gauche"]["ref"] = [
    [1.092020143326, 0.163175911167, -0.925416578398],
    [0.750000000000, 0.000000000000, 0.000000000000],
    [-0.750000000000, 0.000000000000, 0.000000000000],
    [-1.092020143326, -0.163175911167, -0.925416578398],
    [-1.092020143326, -0.719846310393, 0.604022773555],
    [-1.092020143326, 0.883022221559, 0.321393804843],
    [1.092020143326, 0.719846310393, 0.604022773555],
    [1.092020143326, -0.883022221559, 0.321393804843],
]


data["triplet_ethylene"]["ref"] = [
    [0.000000000000, 0.000000000000, -0.705000000000],
    [0.000000000000, 0.000000000000, 0.705000000000],
    [0.000000000000, 0.924372424811, -1.282611998014],
    [0.000000000000, -0.924372424811, -1.282611998014],
    [-0.924372424811, 0.000000000000, 1.282611998014],
    [0.924372424811, 0.000000000000, 1.282611998014],
]


data["allene"]["ref"] = [
    [-2.000000000000, -0.707106781187, 0.707106781187],
    [-2.000000000000, 0.707106781187, -0.707106781187],
    [-1.500000000000, -0.000000000000, 0.000000000000],
    [0.000000000000, -0.000000000000, 0.000000000000],
    [1.500000000000, 0.000000000000, 0.000000000000],
    [2.000000000000, 0.707106781187, 0.707106781187],
    [2.000000000000, -0.707106781187, -0.707106781187],
]


data["ethane_staggered"]["ref"] = [
    [-0.657774170895, 0.990250821217, -0.813797681349],
    [-0.134283313545, 0.737880743551, 0.000000000000],
    [0.134283313545, -0.737880743551, 0.000000000000],
    [0.657774170895, -0.990250821217, -0.813797681349],
    [-0.728988008574, -1.242620898884, 0.000000000000],
    [0.657774170895, -0.990250821217, 0.813797681349],
    [-0.657774170895, 0.990250821217, 0.813797681349],
    [0.728988008574, 1.242620898884, 0.000000000000],
]


data["singlet_ethylene"]["ref"] = [
    [0.000000000000, 0.000000000000, -0.705000000000],
    [0.000000000000, 0.000000000000, 0.705000000000],
    [0.000000000000, 0.924372424811, -1.282611998014],
    [0.000000000000, -0.924372424811, -1.282611998014],
    [0.000000000000, 0.924372424811, 1.282611998014],
    [0.000000000000, -0.924372424811, 1.282611998014],
]


data["ethane_eclipsed"]["ref"] = [
    [0.000000000000, 1.092020143326, -0.939692620786],
    [0.000000000000, 0.750000000000, 0.000000000000],
    [0.000000000000, -0.750000000000, 0.000000000000],
    [0.000000000000, -1.092020143326, -0.939692620786],
    [0.813797681349, -1.092020143326, 0.469846310393],
    [-0.813797681349, -1.092020143326, 0.469846310393],
    [-0.813797681349, 1.092020143326, 0.469846310393],
    [0.813797681349, 1.092020143326, 0.469846310393],
]


data["BH4p"]["ref"] = [
    [0.000000000000, 0.000000000000, -1.000000000000],  # Gh
    [0.000000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, -1.000000000000, 0.000000000000],
    [1.000000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, 1.000000000000, 0.000000000000],
    [-1.000000000000, 0.000000000000, 0.000000000000],
]


data["CH4"]["ref"] = [
    [0.000000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, 0.889981273211, -0.629311793417],
    [0.889981273211, 0.000000000000, 0.629311793417],
    [-0.889981273211, 0.000000000000, 0.629311793417],
    [0.000000000000, -0.889981273211, -0.629311793417],
]


data["SF6"]["ref"] = [
    [0.000000000000, 0.000000000000, -1.800000000000],
    [0.000000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, -1.800000000000, 0.000000000000],
    [1.800000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, 1.800000000000, 0.000000000000],
    [-1.800000000000, 0.000000000000, 0.000000000000],
    [0.000000000000, 0.000000000000, 1.800000000000],
]


data["Ih"]["ref"] = [
    [-0.000000000000, 1.000000000000, 1.618033988750],
    [0.000000000000, -1.000000000000, 1.618033988750],
    [-0.000000000000, 1.000000000000, -1.618033988750],
    [0.000000000000, -1.000000000000, -1.618033988750],
    [1.000000000000, 1.618033988750, 0.000000000000],
    [-1.000000000000, 1.618033988750, 0.000000000000],
    [1.000000000000, -1.618033988750, 0.000000000000],
    [-1.000000000000, -1.618033988750, 0.000000000000],
    [1.618033988750, 0.000000000000, 1.000000000000],
    [1.618033988750, 0.000000000000, -1.000000000000],
    [-1.618033988750, -0.000000000000, 1.000000000000],
    [-1.618033988750, -0.000000000000, -1.000000000000],
]


data["H2_mol"]["ref"] = [[-0.37, 0, 0], [0.37, 0, 0]]

#! Tests to determine full point group symmetry.  Currently, these only matter
#! for the rotational symmetry number in thermodynamic computations.

data["C2H2"]["pg"] = "D_inf_h"
data["C2H2"]["rsn"] = 2
data["C2H2"][
    "mol"
] = """
  {isoA}C 0 0  r1
  C 0 0 -r1
  H 0 0  r2
  H 0 0 -r2
  r1 = 0.65
  r2 = 1.75
"""
data["isoC2H2"]["A"] = 13
data["isoC2H2"]["pg"] = "C_inf_v"
data["isoC2H2"]["rsn"] = 1
data["isoC2H2"]["A"] = 13


data["N2"]["pg"] = "D_inf_h"
data["N2"]["rsn"] = 2
data["N2"][
    "mol"
] = """
  {isoA}N 0.0 0.0 0.0
  N 0.0 0.0 r
  r = 1.1
"""
data["isoN2"]["A"] = 15
data["isoN2"]["pg"] = "C_inf_v"
data["isoN2"]["rsn"] = 1


data["CN"]["pg"] = "C_inf_v"
data["CN"]["rsn"] = 1
data["CN"][
    "mol"
] = """
  0 2
  {isoA}C 0.0 0.0 0.0
  N 0.0 0.0 r
  r = 1.4
"""
data["isoCN"]["A"] = 13
data["isoCN"]["pg"] = "C_inf_v"
data["isoCN"]["rsn"] = 1


data["HCCCl"]["pg"] = "C_inf_v"
data["HCCCl"]["rsn"] = 1
data["HCCCl"][
    "mol"
] = """
  H  0 0 -1.0
  {isoA}C  0 0  0.0
  C  0 0  1.1
  Cl 0 0  2.6
"""
data["isoHCCCl"]["A"] = 14
data["isoHCCCl"]["pg"] = "C_inf_v"
data["isoHCCCl"]["rsn"] = 1


data["CHFClBr"]["pg"] = "C1"
data["CHFClBr"]["rsn"] = 1
data["CHFClBr"][
    "mol"
] = """
  {isoA}H
  C  1 1.0
  F  2 1.0 1 105.0
  Cl 2 1.0 1 105.0 3  120.0
  Br 2 1.0 1 105.0 3 -120.0
"""
data["isoCHFClBr"]["A"] = 2
data["isoCHFClBr"]["pg"] = "C1"
data["isoCHFClBr"]["rsn"] = 1


data["CH2ClBr"]["pg"] = "Cs"
data["CH2ClBr"]["rsn"] = 1
data["CH2ClBr"][
    "mol"
] = """
  Cl
  C  1 1.0
  Br 2 1.0 1 105.0
  H  2 1.0 1 105.0 3  120.0
  {isoA}H  2 1.0 1 105.0 3 -120.0
"""
data["isoCH2ClBr"]["A"] = 2
data["isoCH2ClBr"]["pg"] = "C1"
data["isoCH2ClBr"]["rsn"] = 1


data["HOCl"]["pg"] = "Cs"
data["HOCl"]["rsn"] = 1
data["HOCl"][
    "mol"
] = """
  {isoA}H
  O 1 1.0
  Cl 2 1.7 1 110.0
"""
data["isoHOCl"]["A"] = 2
data["isoHOCl"]["pg"] = "Cs"
data["isoHOCl"]["rsn"] = 1


data["C4H4Cl2F2"]["au"] = True
data["C4H4Cl2F2"]["pg"] = "Ci"
data["C4H4Cl2F2"]["rsn"] = 1
data["C4H4Cl2F2"][
    "mol"
] = """
  units    bohr
  C     0.432781050498     1.898774028282     0.810337938486
  C    -1.658744642774     0.805191018766    -0.984829058337
  {isoA}C     1.658744642774    -0.805191018766     0.984829058337
  C    -0.432781050498    -1.898774028282    -0.810337938486
  H    -0.317971784026     2.532165941971     2.640915161238
  H    -1.615729990528     1.614062700629    -2.881498569657
  H     1.615729990528    -1.614062700629     2.881498569657
  H     0.317971784026    -2.532165941971    -2.640915161238
  Cl   -4.852178875691     1.024620478757     0.190249941464
  Cl    4.852178875691    -1.024620478757    -0.190249941464
  F    -1.913713787211    -3.739567959534     0.258534542158
  F     1.913713787211     3.739567959534    -0.258534542158
"""
data["isoC4H4Cl2F2"]["A"] = 13
data["isoC4H4Cl2F2"]["pg"] = "C1"
data["isoC4H4Cl2F2"]["rsn"] = 1


data["HOOH_dimer"]["pg"] = "Ci"
data["HOOH_dimer"]["rsn"] = 1
data["HOOH_dimer"][
    "mol"
] = """
  H   0.9911262285  -1.7979226333   0.1465182515
  O   2.7691093095  -1.3485218649  -0.0071557684
  O   2.5178030311   1.3808374923  -0.1154058014
  H   3.2883200453   1.8308595095   1.4757706825
  H  -3.2883200453  -1.8308595095  -1.4757706825
  O  -2.5178030311  -1.3808374923   0.1154058014
  O  -2.7691093095   1.3485218649   0.0071557684
  {isoA}H  -0.9911262285   1.7979226333  -0.1465182515
"""
data["isoHOOH_dimer"]["A"] = 2
data["isoHOOH_dimer"]["pg"] = "C1"
data["isoHOOH_dimer"]["rsn"] = 1


data["HOOH"]["pg"] = "C2"
data["HOOH"]["rsn"] = 2
data["HOOH"][
    "mol"
] = """
  H
  O 1 1.0
  {isoA}O 2 1.5 1 110.0
  H 3 1.0 2 110.0 1 60.0
"""
data["isoHOOH"]["A"] = 18
data["isoHOOH"]["pg"] = "C1"
data["isoHOOH"]["rsn"] = 1


data["NOHOHOH"]["pg"] = "C3"
data["NOHOHOH"]["rsn"] = 3
data["NOHOHOH"][
    "mol"
] = """
  X
  N 1 1.0
  {isoA}O 2 1.5 1 100.0
  O 2 1.5 1 100.0  3  120.0
  O 2 1.5 1 100.0  3 -120.0
  H 3 1.0 2 110.0 4 0.0
  H 4 1.0 2 110.0 5 0.0
  H 5 1.0 2 110.0 3 0.0
"""
data["isoNOHOHOH"]["A"] = 18
data["isoNOHOHOH"]["pg"] = "Cs"
data["isoNOHOHOH"]["rsn"] = 1


data["H2O"]["pg"] = "C2v"
data["H2O"]["rsn"] = 2
data["H2O"][
    "mol"
] = """
  H
  O 1 1.0
  {isoA}H 2 1.0 1 109.5
"""
data["isoH2O"]["A"] = 2
data["isoH2O"]["pg"] = "Cs"
data["isoH2O"]["rsn"] = 1


data["CH2F2"]["au"] = True
data["CH2F2"]["pg"] = "C2v"
data["CH2F2"]["rsn"] = 2
data["CH2F2"][
    "mol"
] = """
  units au
  C     0.0000000000  -0.0000000000   1.0890958457
  F     0.0000000000  -2.1223155812  -0.4598161475
  F    -0.0000000000   2.1223155812  -0.4598161475
  {isoA}H     1.7084139850   0.0000000000   2.1841068002
  H    -1.7084139850  -0.0000000000   2.1841068002
"""
data["isoCH2F2"]["A"] = 3
data["isoCH2F2"]["pg"] = "Cs"
data["isoCH2F2"]["rsn"] = 1


data["NH3"]["pg"] = "C3v"
data["NH3"]["rsn"] = 3
data["NH3"][
    "mol"
] = """
  X
  N 1 1.0
  H 2 rNH 1 aXNH
  {isoA}H 2 rNH 1 aXNH 3 120.0
  H 2 rNH 1 aXNH 4 120.0

  rNH = 0.95
  aXNH = 115.0
"""
data["isoNH3"]["A"] = 3
data["isoNH3"]["pg"] = "Cs"
data["isoNH3"]["rsn"] = 1


data["BrF5"]["pg"] = "C4v"
data["BrF5"]["rsn"] = 4
data["BrF5"][
    "mol"
] = """
 F
 Br 1 r
 F  2 r 1 90.0
 {isoA}F  2 r 3 90.0 1  90.0
 F  2 r 3 90.0 1 -90.0
 F  2 r 1 90.0 3 180.0
 r = 1.7
"""
data["isoBrF5"]["A"] = 18
data["isoBrF5"]["pg"] = "Cs"
data["isoBrF5"]["rsn"] = 1


data["N2H2"]["pg"] = "C2h"
data["N2H2"]["rsn"] = 2
data["N2H2"][
    "mol"
] = """
  N
  N 1 rNN
  H 1 rNH 2 aHNN
  H 2 rNH 1 aHNN 3 180.0
  rNH  = 1.0
  rNN  = 1.4
  aHNN = 140.0
"""


data["NOHOHOH"]["pg"] = "C3h"
data["NOHOHOH"]["rsn"] = 3
data["NOHOHOH"][
    "mol"
] = """
  X
  N 1 1.0
  O 2 1.5 1 90.0
  O 2 1.5 1 90.0  3  120.0
  O 2 1.5 1 90.0  3 -120.0
  {isoA}H 3 1.0 2 110.0 4 0.0
  H 4 1.0 2 110.0 5 0.0
  H 5 1.0 2 110.0 3 0.0
"""
data["isoNOHOHOH"]["A"] = 2
data["isoNOHOHOH"]["pg"] = "Cs"
data["isoNOHOHOH"]["rsn"] = 1


# 1,3,5,7-tetrafluorocyclooctatetraene
data["TFCOT"]["pg"] = "S4"
data["TFCOT"]["rsn"] = 2
data["TFCOT"][
    "mol"
] = """
  C       -1.618188     -0.437140     -0.409373
  C       -1.394411      0.896360     -0.429596
  C       -0.896360     -1.394411      0.429596
  C       -0.437140      1.618188      0.409373
  C        0.437140     -1.618188      0.409373
  C        0.896360      1.394411      0.429596
  C        1.394411     -0.896360     -0.429596
  C        1.618188      0.437140     -0.409373
  F        2.147277     -1.690111     -1.235043
  F        1.690111      2.147277      1.235043
  F       -2.147277      1.690111     -1.235043
  F       -1.690111     -2.147277      1.235043
  H        0.878010     -2.418132      1.029595
  H       -2.418132     -0.878010     -1.029595
  H       -0.878010      2.418132      1.029595
  H        2.418132      0.878010     -1.029595
"""


data["Li_H2O_4_p"]["pg"] = "S4"
data["Li_H2O_4_p"]["rsn"] = 2
data["Li_H2O_4_p"][
    "mol"
] = """
   1 1
   X
   Li 1 1.0
   X 2 1.0 1 90.0
   X 2 1.0 3 90.0 1 180.0
   O 2 oli 1 olix 3 -90.0
   O 2 oli 1 olix 3 90.0
   O 2 oli 4 olix 3 0.0
   O 2 oli 4 olix 3 180.0
   H 5 oh1 2 lioh1 1 xlioh1
   H 5 oh2 2 lioh2 1 xlioh2
   H 6 oh1 2 lioh1 1 xlioh1
   H 6 oh2 2 lioh2 1 xlioh2
   H 7 oh1 2 lioh1 4 -xlioh1
   H 7 oh2 2 lioh2 4 -xlioh2
   H 8 oh1 2 lioh1 4 -xlioh1
   H 8 oh2 2 lioh2 4 -xlioh2
   olix=52.0
   oli=1.9
   oh1=0.952
   oh2=0.950
   lioh1=125.4
   lioh2=124.8
   xlioh1=-40.0
   xlioh2=135.0
"""


data["ethylene_cation"]["pg"] = "D2"
data["ethylene_cation"]["rsn"] = 4
data["ethylene_cation"][
    "mol"
] = """
  C1
  C2 C1 rCC
  H1 C1 rCH C2 aHCC
  H2 C1 rCH C2 aHCC H1 180.0
  H3 C2 rCH C1 aHCC H1 D
  H4 C2 rCH C1 aHCC H3 180.0
  rCC  = 1.41
  rCH  = 1.09
  aHCC = 122.0
  D    = 45.0
"""


data["ethane_gauche"]["pg"] = "D3"
data["ethane_gauche"]["rsn"] = 6
data["ethane_gauche"][
    "mol"
] = """
  H
  C 1 1.0
  C 2 1.5 1 110.0
  H 3 1.0 2 110.0 1   20.0
  H 3 1.0 2 110.0 1  140.0
  H 3 1.0 2 110.0 1 -100.0
  H 2 1.0 3 110.0 1  120.0
  H 2 1.0 3 110.0 1 -120.0
"""


data["triplet_ethylene"]["pg"] = "D2d"
data["triplet_ethylene"]["rsn"] = 4
data["triplet_ethylene"][
    "mol"
] = """
  0 3
  C1
  C2 C1 rCC
  H1 C1 rCH C2 aHCC
  H2 C1 rCH C2 aHCC H1 180.0
  H3 C2 rCH C1 aHCC H1 D
  H4 C2 rCH C1 aHCC H3 180.0
  rCC  = 1.41
  rCH  = 1.09
  aHCC = 122.0
  D    = 90.0
"""


data["allene"]["pg"] = "D2d"
data["allene"]["rsn"] = 4
data["allene"][
    "mol"
] = """
  {isoA}H -2.0  0.0  1.0
  H -2.0  0.0 -1.0
  C -1.5  0.0  0.0
  C  0.0  0.0  0.0
  C  1.5  0.0  0.0
  H  2.0  1.0  0.0
  {isoA}H  2.0 -1.0  0.0
"""
data["isoallene"]["A"] = 2
data["isoallene"]["pg"] = "C2"
data["isoallene"]["rsn"] = 2


data["ethane_staggered"]["pg"] = "D3d"
data["ethane_staggered"]["rsn"] = 6
data["ethane_staggered"][
    "mol"
] = """
  H
  C 1 1.0
  C 2 1.5 1 110.0
  H 3 1.0 2 110.0 1   60.0
  H 3 1.0 2 110.0 1  -60.0
  H 3 1.0 2 110.0 1  180.0
  H 2 1.0 3 110.0 1  120.0
  H 2 1.0 3 110.0 1 -120.0
"""


data["singlet_ethylene"]["pg"] = "D2h"
data["singlet_ethylene"]["rsn"] = 4
data["singlet_ethylene"][
    "mol"
] = """
    C1
    C2 C1 rCC
    {isoA}H1 C1 rCH C2 aHCC
    H2 C1 rCH C2 aHCC H1 180.0
    {isoA}H3 C2 rCH C1 aHCC H1 D
    H4 C2 rCH C1 aHCC H3 180.0
    rCC  = 1.41
    rCH  = 1.09
    aHCC = 122.0
    D    = 0.0
"""
data["isosinglet_ethylene"]["A"] = 2
data["isosinglet_ethylene"]["pg"] = "C2v"
data["isosinglet_ethylene"]["rsn"] = 2


data["ethane_eclipsed"]["pg"] = "D3h"
data["ethane_eclipsed"]["rsn"] = 6
data["ethane_eclipsed"][
    "mol"
] = """
  H
  C 1 1.0
  C 2 1.5 1 110.0
  H 3 1.0 2 110.0 1   00.0
  H 3 1.0 2 110.0 1  120.0
  H 3 1.0 2 110.0 1 -120.0
  H 2 1.0 3 110.0 1  120.0
  H 2 1.0 3 110.0 1 -120.0
"""


data["BH4p"]["pg"] = "D4h"
data["BH4p"]["rsn"] = 8
data["BH4p"][
    "mol"
] = """
 1 1
 X
 B 1 1.0
 {isoA}H 2 1.0 1 90.0
 H 2 1.0 1 90.0 3  90.0
 {isoA}H 2 1.0 1 90.0 3 180.0
 H 2 1.0 1 90.0 3 -90.0
"""
data["isoBH4p"]["A"] = 2
data["isoBH4p"]["pg"] = "D2h"
data["isoBH4p"]["rsn"] = 4


data["CH4"]["pg"] = "Td"
data["CH4"]["rsn"] = 12
data["CH4"][
    "mol"
] = """
   C
   {isoA}H 1 r
   H 1 r 2 TDA
   H 1 r 2 TDA 3 120
   H 1 r 2 TDA 4 120
   r = 1.09
"""
data["isoCH4"]["A"] = 2
data["isoCH4"]["pg"] = "C3v"
data["isoCH4"]["rsn"] = 3


data["SF6"]["pg"] = "Oh"
data["SF6"]["rsn"] = 24
data["SF6"][
    "mol"
] = """
  {isoA}F
  S 1 r
  {isoA}F 2 r 1 90.0
  F 2 r 1 90.0 3  90.0
  {isoA}F 2 r 1 90.0 3 180.0
  F 2 r 1 90.0 3 -90.0
  {isoA}F 2 r 5 90.0 1 180.0
  r = 1.8
"""
data["isoSF6"]["A"] = 18
data["isoSF6"]["pg"] = "D4h"
data["isoSF6"]["rsn"] = 8


data["Ih"]["au"] = True
data["Ih"]["pg"] = "Ih"
data["Ih"]["rsn"] = 60
data["Ih"][
    "mol"
] = """
  unit = au
  0 1
  H   0   1   x
  {isoA}H   0  -1   x
  H   0   1  -x
  H   0  -1  -x
  {isoA}H   1   x   0
  {isoA}H  -1   x   0
  H   1  -x   0
  H  -1  -x   0
  {isoA}H   x   0   1
  H   x   0  -1
  {isoA}H  -x   0   1
  H  -x   0  -1
  x = 1.618033988749894848
"""
data["isoIh"]["A"] = 2
data["isoIh"]["pg"] = "C5v"
data["isoIh"]["rsn"] = 5


data["H2_mol"]["pg"] = "D_inf_h"
data["H2_mol"]["rsn"] = 2
data["H2_mol"]["mol"] = """H 0 0 0\nH 0.74 0 0"""


@pytest.mark.parametrize("subject", data.keys())
def test_molsymm_alone(subject):
    pg = data[subject]["pg"]
    sigma = data[subject]["rsn"]

    if subject.startswith("iso"):
        isbohr = data[subject[3:]].get("au", False)
        molstr = data[subject[3:]]["mol"].format(isoA=data[subject]["A"])
        refgeomang = None
    else:
        isbohr = data[subject].get("au", False)
        molstr = data[subject]["mol"].format(isoA="")
        refgeomang = data[subject]["ref"]

    symmol = qcdb.Molecule(molstr)
    symmol.update_geometry()
    # symmol.axis_representation()

    assert compare(pg, symmol.get_full_point_group(), pg + " point group: " + subject)
    assert compare(sigma, symmol.rotational_symmetry_number(), pg + " sigma")

    if isbohr:
        geom_now = symmol.full_geometry()
    else:
        geom_now = qcdb.util.vecutil.mscale(symmol.full_geometry(), qcel.constants.bohr2angstroms)

    if refgeomang:
        assert compare_values(refgeomang, geom_now, pg + " orientation", atol=1.0e-6)


@pytest.mark.parametrize("subject", [subject for subject in data.keys() if not subject.startswith("iso")])
@pytest.mark.parametrize(
    "qcprog",
    [
        pytest.param("cfour", marks=using("cfour")),
        pytest.param("gamess", marks=using("gamess")),
        pytest.param("nwchem", marks=using("nwchem")),
        pytest.param("psi4", marks=using("psi4")),
    ],
)
def test_molsymm_calc(subject, qcprog):

    from .data_mints5 import mints5

    ref_block = mints5[subject]
    expected = data[subject]
    abbr = {"cfour": "c4", "gamess": "gms", "nwchem": "nwc", "psi4": "p4"}

    symmol = qcdb.Molecule(data[subject]["mol"].format(isoA=""))  # since no iso subjects yet. later: expected["A"]
    symmol.update_geometry()

    assert compare(expected["pg"], symmol.get_full_point_group(), "point group")
    assert compare(expected["rsn"], symmol.rotational_symmetry_number(), "sigma")
    units = "bohr" if expected.get("au", False) else "angstrom"
    refau = np.array(data[subject]["ref"]) * qcel.constants.conversion_factor(units, "bohr")
    assert compare_values(refau, symmol.full_geometry(), atol=1.0e-6)

    if subject == "Ih":
        if qcprog in ["cfour", "gamess"]:
            pytest.xfail("weird hang, probably atom permutations")
        elif qcprog == "nwchem":
            pytest.xfail("weird gradient rotation")
    elif subject == "SF6" and qcprog == "nwchem":
        pytest.xfail(
            "TODO new bug 'autosym bug : too many atoms:Received an Error in Communication' after fixed basis set for ghosts"
        )

    extra_keywords = {
        "cfour": {
            "HOOH_dimer": {"cfour_scf_damping": 900},
            "TFCOT": {"cfour_scf_damping": 800},
        },
        "gamess": {
            "HOOH_dimer": {"gamess_contrl__maxit": 60, "gamess_scf__conv": 1.0e-7},
        },
        "nwchem": {
            "CN": {"nwchem_scf__thresh": 1.0e-8},
            "HOOH_dimer": {"nwchem_scf__nr": 0.0, "nwchem_scf__maxiter": 100, "nwchem_scf__thresh": 1.0e-6},
        },
        "psi4": {
            "HOOH_dimer": {"psi4_soscf": True},
            "Ih": {"psi4_guess": "gwh"},
        },
    }
    qcdb.set_options(
        {
            "scf_type": "pk",
            "e_convergence": 9,
            "d_convergence": 9,
            "reference": "r{}hf".format("" if symmol.multiplicity() == 1 else "o"),
        }
    )
    qcdb.set_options(extra_keywords[qcprog].get(subject, {}))

    model = f"{abbr[qcprog]}-hf/cc-pvdz"
    grad, wfn = qcdb.gradient(model, return_wfn=True, molecule=symmol, local_options={"ncores": 2})
    pprint.pprint(wfn, width=200)

    with np.printoptions(precision=3, suppress=True):
        # DEBUG
        pprint.pprint(ref_block)
        print("refau")
        print(refau)
        print("full_geometry")
        print(symmol.full_geometry(np_out=True))
        print("wfn mol")
        print(wfn["molecule"]["geometry"])
        print("frad")
        print(grad)
        dgrad = np.array2string(grad, separator=", ")
        print(dgrad)

    for i, item in enumerate(
        [
            symmol.nuclear_repulsion_energy(),
            wfn["qcvars"]["NUCLEAR REPULSION ENERGY"].data,
            wfn["extras"]["qcvars"]["NUCLEAR REPULSION ENERGY"],
            wfn["properties"]["nuclear_repulsion_energy"],
        ]
    ):
        assert compare_values(
            ref_block["nre"], item, atol=1.0e-6, label="NRE"
        ), f"NRE {i}: {item} != (ref) {ref_block['nre']}"

    for i, item in enumerate(
        [
            wfn["qcvars"]["HF TOTAL ENERGY"].data,
            wfn["qcvars"]["CURRENT ENERGY"].data,
            wfn["extras"]["qcvars"]["HF TOTAL ENERGY"],
            wfn["extras"]["qcvars"]["CURRENT ENERGY"],
            wfn["properties"]["return_energy"],
            wfn["properties"]["scf_total_energy"],
        ]
    ):
        assert compare_values(
            ref_block["hf_total_energy"], item, atol=1.0e-6, label="ene"
        ), f"Energy {i}: {item} != (ref) {ref_block['hf_total_energy']}"

    assert compare(ref_block["elez"], wfn["molecule"]["atomic_numbers"], label="elez"), "elez"
    assert compare(ref_block["elea"], wfn["molecule"]["mass_numbers"], label="elea"), "elea"
    assert compare(False, wfn["molecule"]["fix_orientation"], label="frame"), "frame"

    for i, item in enumerate(
        [
            symmol.geometry(),
            wfn["molecule"]["geometry"],
            # qcdb.driver.driver_helpers.get_active_molecule().geometry(),
        ]
    ):
        assert compare_values(
            np.array(ref_block["geom"]).reshape((-1, 3)), item, atol=1.0e-6, label="geom"
        ), f"Geometry {i}"

    for i, item in enumerate(
        [
            grad,
            wfn["qcvars"]["HF TOTAL GRADIENT"].data,
            wfn["qcvars"]["CURRENT GRADIENT"].data,
            wfn["extras"]["qcvars"]["HF TOTAL GRADIENT"],
            wfn["extras"]["qcvars"]["CURRENT GRADIENT"],
            # wfn["properties"]["return_gradient"],
            np.array(wfn["properties"]["scf_total_gradient"]).reshape((-1, 3)),  # TODO shouldn't need to reshape
        ]
    ):
        assert compare_values(
            np.array(ref_block["hf_total_gradient"]).reshape((-1, 3)), item, atol=1.0e-6, label="grad"
        ), f"Gradient {i}"

    # assert compare_values(np.array(ref_block["hf_total_dipole"]), dipole, atol=1.e-6, label="dipole")
    # assert 0


# Ih:      options.add_str("GUESS", "AUTO", "AUTO CORE GWH SAD SADNO SAP HUCKEL READ");
#           '             Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u \n'
#           '    DOCC [     2,    0,    0,    1,    0,    1,    1,    1 ]\n' '  @RHF Final Energy:    -5.45327901072954\n' sad & huckel
#           '    DOCC [     1,    1,    1,    0,    0,    1,    1,    1 ]\n' '  @RHF Final Energy:    -5.51322075413067\n' sap
#           '    DOCC [     2,    1,    0,    0,    0,    1,    1,    1 ]\n' '  @RHF Final Energy:    -5.53971699010802\n' core & sadno
#           '    DOCC [     3,    0,    0,    0,    0,    1,    1,    1 ]\n' '  @RHF Final Energy:    -5.54050178487978\n' gwh
