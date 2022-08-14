# TomoNV example, Python version, case 1.  information for specific orientation

from tomoNV_Python import *

InitPARAMS(
  bUseRAY = False,        #use multi-thrading for fast calculation. not recommeded here.
  bUseNUMBA = True,       #use NUMBA partially for fast calculation
  bUseSlotPairing = True, #for debugging. 
  bUseExplicitSS = False) #Set bUseExplicitSS=True to see the SS pixels visually.

TestFiles = ["TomoNV_Data\\(6)Bunny_5k.stl"]

import os.path
time0 = StartTimer()
for file in TestFiles:
  g_input_mesh_filename = file
  if( os.path.isfile(g_input_mesh_filename) ):
    tomo1 = tomoNV( g_input_mesh_filename, bVerbose=True)

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
    tomo1.yaw   = toRadian(0.) #initial orientation. change these when you want to rotate the input mesh manually
    tomo1.pitch = toRadian(0.) #initial orientation. change these when you want to rotate the input mesh manually
    tomo1.roll  = toRadian(0.) #initial orientation. change these when you want to rotate the input mesh manually

    tomo1.Pixelize()
    tomo1.Calculate()
    tomo1.Print()# print final volume/mass information.

    # Plot2DTomo (tomo1) #show 2D tomographs
    Plot3DPixels( tomo1, tomo1.FigureTitle()) #show 3D pixles
    #tomo1.Print_tabbed(tomo1)#print mass info as tab-separated string, for debugging
    #PrintSlotInfo( tomo1, X=29,Y=29 ) #print pxls in (X,Y) slot. for debugging
  else:
    print(Fore.RED , '[ERROR] File does not exists: ' , g_input_mesh_filename ,  Style.RESET_ALL)
EndTimer( time0)
