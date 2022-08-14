# TomoNV example, Python version, Case 2. optimal orientation searching

from tomoNV_Python import *

InitPARAMS(
  bUseRAY = True,         #use multi-thrading for fast calculation. "True" is recommended here.
  bUseNUMBA = True,       #use NUMBA partially for fast calculation
  bUseSlotPairing = True, #for debugging. 
  bUseExplicitSS = False) ##Set bUseExplicitSS=True to see the SS pixels visually.

g_input_mesh_filename = 'TomoNV_Data\\(6)Bunny_5k.stl'

import os.path
import sys
if( not os.path.isfile(g_input_mesh_filename) ):
  print(Fore.RED + '[ERROR] File Does not Exists: ' + g_input_mesh_filename +  Style.RESET_ALL)
  sys.exit(0)

#file loading
tomo1 = tomoNV( g_input_mesh_filename, bVerbose=False)

'''
#default g-code parameters
tomo1.wall_thickness = 0.8 # [mm]
tomo1.PLA_density    = 0.00121 # density of PLA filament, [g/mm^3]
tomo1.Fclad = 1.0 # fill ratio of cladding, always 1.0
tomo1.Fcore   = 0.15 # fill ratio of core, (0~1.0)
tomo1.Fss     = 0.2 # fill ratio of support structure, (0~1.0)
tomo1.dVoxel = 1.0
tomo1.nVoxel = 256
'''

#user options
tomo1.Css   = 1.2 # correction constant for filament dilation effect. 
tomo1.theta_c = toRadian(60.) #support structure critical angle
theta_YP = 45 #angle interval. in degree, in integer. The smaller value takes the more time.

#search optimal
nYPR_Intervals= int(360 / theta_YP) +1
yaw_range   = np.linspace(toRadian(0), toRadian(360), num=nYPR_Intervals, endpoint=True, dtype=np.float32)
pitch_range = np.linspace(toRadian(0), toRadian(360), num=nYPR_Intervals, endpoint=True, dtype=np.float32)
roll  = np.zeros(1)

time0 = StartTimer()
SearchYP(  tomo1, yaw_range, pitch_range, bVerbose=True, bShowGraph=True) # draw contour graphs, optimal and worst orientations.
EndTimer( time0)

i=0