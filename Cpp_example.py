# TomoNV example, C++ Dll version
from  tomoNV_Cpp import *

#-----------------------------------------------
#(1) specify (filename, initial orientation) and angle interval.
#========================================================================================================================
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(1)cone50_63k.obj',     0, 0, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(1)cone50_63k.obj',    10, 332, 0) #1st worst

# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(2)sphere90_39k.obj',    0, 0, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(2)sphere90_39k.obj',  180, 270, 0) #1st worst

# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(3)Bunny_69k_2x.obj', 231,  50, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(3)Bunny_69k_2x.obj', 243, 357, 0) #1st worst

(g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(4)Bunny_69k.stl', 247, 46, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(4)Bunny_69k.stl', 243, 346, 0) #1st worst

# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(5)Bunny_69k_0.5x.stl',  72, 130, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(5)Bunny_69k_0.5x.stl',  51, 200, 0) #1st worst

# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(6)Bunny_5k.stl',  239,  52, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(6)Bunny_5k.stl',  63, 182, 0) #1st worst

# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(7)Bunny_1k.obj', 284,  48, 0) #1st optimal
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(7)Bunny_1k.obj', 241, 348, 0) #1st worst

# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(8)manikin.obj',  139, 277, 0) #1st worst
# (g_input_mesh_filename, Yaw, Pitch, Roll) = ('TomoNV_Data\\(8)manikin.obj',  268,  20, 0) #1st optimal

theta_YP = 10 #angle iterval. should be "360 / integer". if 0, see the specified initial orientation only.
#========================================================================================================================
  
if(theta_YP==0):
  # case 1). seeing mass infomatin of a specific orientation.  
  nYPR_Intervals=1
  yaw_range   = np.ones(nYPR_Intervals) * toRadian(Yaw)
  pitch_range = np.ones(nYPR_Intervals) * toRadian(Pitch)
  roll_range  = np.ones(nYPR_Intervals) * toRadian(Roll)
elif(theta_YP> 0):
  # case 2). searching optimal orientation
  nYPR_Intervals= int(360 / theta_YP) +1
  yaw_range   = np.linspace(toRadian(0), toRadian(360), num=nYPR_Intervals, endpoint=True, dtype=np.float32)#change angle interval manually for adaptive search.
  pitch_range = np.linspace(toRadian(0), toRadian(360), num=nYPR_Intervals, endpoint=True, dtype=np.float32)
  roll_range  = np.zeros(1) #for generality. roll direction is not needed in our modeling.
else:
  import sys
  sys.exit(0)

tomoNV_Cpp1 = tomoNV_Cpp(g_input_mesh_filename, nYPR_Intervals, yaw_range, pitch_range, roll_range )# load mesh file and Cpp Dll
#-----------------------------------------------
# (2) input g-code conditions and run C++ DLL.

'''
#default g-code parameters
tomoNV_Cpp1.wall_thickness = 0.8 # [mm]
tomoNV_Cpp1.PLA_density    = 0.00121 # density of PLA filament, [g/mm^3]
tomoNV_Cpp1.Fclad = 1.0 # fill ratio of cladding, always 1.0
tomoNV_Cpp1.Fcore   = 0.15 # fill ratio of core, (0~1.0)
tomoNV_Cpp1.Fss     = 0.2 # fill ratio of support structure, (0~1.0)
tomoNV_Cpp1.dVoxel = 1.0
tomoNV_Cpp1.nVoxel = 256
'''

#user options
tomoNV_Cpp1.Css   = 1.2 # correction constant for filament dilation effect. 
tomoNV_Cpp1.theta_c = toRadian(60.)#support structure critical angle
tomoNV_Cpp1.bUseExplicitSS = True #Set as True to see the SS pixels visually.

time0 = StartTimer()
tomoNV_Cpp1.Run(cpp_function_name = 'TomoNV_TMPxl') #call C++ 
EndTimer( time0)
#-----------------------------------------------
# (3) Rendering

Plot3DPixels(tomoNV_Cpp1) #show pixels of the 1st optimal

if tomoNV_Cpp1.nYPR_Intervals >= 5:
  (optimal_YPRs,worst_YPRs) = findOptimals(tomoNV_Cpp1.YPR, tomoNV_Cpp1.Mss3D, g_nOptimalsToDisplay) #find optimal and worst orientation
  Plot3D(tomoNV_Cpp1.mesh0, yaw_range, pitch_range, tomoNV_Cpp1.Mss3D, optimal_YPRs, worst_YPRs) # plot the contour graph like Fig. 20
#-----------------------------------------------
