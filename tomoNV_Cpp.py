# interfaces to call C++ Dll of TomoNV algorithm
from  tomoNV_io import *

class tomoNV_Cpp:
  FileName = ""; mesh0 = [];vtx = []; tri=[]; tri_nrm=[]; xx=[]; yy=[];zz=[];
  YPR=[];DLL_Info=[];CppFunction = [];Cdll_opt_id=[];CppDLL=[];
  Mo3D=[];Mss3D=[];Mtotal3D=[];AABB2D=(0,0,0,0);chull_vtx=[];tri_area=[]
  vm_info=[];Vtc_mm = []

  wall_thickness = 1.0 # [mm]
  PLA_density    = 0.00121 # density of PLA filament, [g/mm^3]
  Fclad = 1.0 # fill ratio of cladding, always 1.0
  Fcore = 0.15 # fill ratio of core, (0~1.0)
  Fss   = 0.2 # fill ratio of support structure, (0~1.0)
  Css   = 1. # correction constant for Mss. obsolete. replaced by getWeightZSum[]
  Fclad = 1.
  bVerbose = True
  bUseExplicitSS = False
  dVoxel = 1 #default: 1mm voxel
  nVoxel = 256 #default: 256 x 256 x 256 voxels
  theta_c = toRadian(60.)
  nYPR_Intervals = 0 #number of intervals. = int(360 / theta_YP) +1
    

  def __init__(self, filename, nYPR_Intervals, yaw_range, pitch_range, roll_range):
    self.FileName = filename
    (self.mesh0, self.vtx, self.tri_nrm, self.tri, self.tri_area, self.tri_center, self.chull_vtx) = LoadInputMesh( filename)
    (self.xx,self.yy,self.zz)   = np.meshgrid( yaw_range, pitch_range, roll_range)
    self.YPR = np.column_stack([self.xx.ravel(), self.yy.ravel(), self.zz.ravel()]).astype(np.float32) #https://rfriend.tistory.com/352
    self.nYPR_Intervals = nYPR_Intervals

    (mesh1_min, mesh1_max,_) = getBoundary( self.vtx)
    (x1,y1,z1) = list(map(int, mesh1_max))
    (x0,y0,z0) = list(map(int, mesh1_min))
    max_dimension = max(max( x1-x0, y1-y0), z1-z0)
    self.AABB2D = (0, 0, max_dimension, max_dimension)#for plot2D()
     
    print('YPR=',self.YPR, self.YPR.shape)  

  def Run(self, cpp_function_name):
    global g_mesh0_surface_area,g_CppDLLFileName
    g_mesh0_surface_area = o3d.geometry.TriangleMesh.get_surface_area( self.mesh0)
    self.DLL_Info = np.array([ self.bVerbose + g_fMARGIN, self.bUseExplicitSS+ g_fMARGIN, self.dVoxel, self.nVoxel + g_fMARGIN, self.YPR.shape[0] + g_fMARGIN, 
      self.theta_c,  self.vtx.shape[0] + g_fMARGIN,  self.tri.shape[0] + g_fMARGIN,  g_mesh0_surface_area,  self.wall_thickness,  
      self.PLA_density,  self.Fcore, self.Fclad,   self.Fss, self.Css,
      self.chull_vtx.shape[0] + g_fMARGIN ]).astype( np.float32) 

    self.CppDLL =  ct.WinDLL(g_CppDLLFileName)    
    self.CppFunction = getattr(self.CppDLL, cpp_function_name)
    self.CppFunction.argtypes   = (  Cptr1d, Cptr1d, Cptr1d, Cptr1d, Cptr1d, Cptr1iL)
    self.CppFunction.restype    = ct.c_int32
    self.Cdll_opt_id = self.CppFunction( 
        np_to_Cptr1d(self.DLL_Info), 
        np_to_Cptr1d(self.YPR),
        np_to_Cptr1d(self.chull_vtx),#chull_vtx is not necessary. 
        np_to_Cptr1d(self.vtx), 
        np_to_Cptr1d(self.tri_nrm), 
        np_to_Cptr1iL(self.tri)  )

    # (5) Retreive data from C++ dll
    self.CppDLL.getMo.argtypes  = (); self.CppDLL.getMo.restype  = Cptr1d
    self.CppDLL.getMss.argtypes = (); self.CppDLL.getMss.restype = Cptr1d

    self.Mo3D  = np.array(Cptr1d_to_np(self.CppDLL.getMo(),  self.YPR.shape[0])).astype(np.float32)
    self.Mss3D = np.array(Cptr1d_to_np(self.CppDLL.getMss(), self.YPR.shape[0])).astype(np.float32)
    self.Mtotal3D = self.Mo3D + self.Mss3D

    self.CppDLL.getVolMassInfo.argtypes = (); self.CppDLL.getVolMassInfo.restype = Cptr1d
    self.vm_info = np.array(Cptr1d_to_np(self.CppDLL.getVolMassInfo(),  13)).astype(np.float32)

    if(self.bVerbose):
      print( Fore.BLUE, 'Mo3D=    ', Style.RESET_ALL, self.Mo3D.reshape(     self.nYPR_Intervals,self.nYPR_Intervals) ) 
      print( Fore.BLUE, 'Mss3D=   ', Style.RESET_ALL, self.Mss3D.reshape(    self.nYPR_Intervals,self.nYPR_Intervals) ) 
      print( Fore.BLUE, 'Mtotal3D=', Style.RESET_ALL, self.Mtotal3D.reshape( self.nYPR_Intervals,self.nYPR_Intervals) ) 
      # np.savetxt('Mo.txt',  Mo3D.reshape( self.nYPR_Intervals,self.nYPR_Intervals), delimiter='\t')
      # np.savetxt('Mss.txt', Mss3D.reshape(self.nYPR_Intervals,self.nYPR_Intervals), delimiter='\t')

    # (6) rendering
    import this
    self.CppDLL.getnData2i.argtypes = ( ct.c_short,); self.CppDLL.getnData2i.restype = ct.c_int32
    self.CppDLL.getpData2i.argtypes = ( ct.c_short,); self.CppDLL.getpData2i.restype = Cptr1i

    for p_type, p_name in zip ( g_PixelEnums, g_PixelVarNames ):
      p_enum = getattr( enumPixelType, p_type).value 
      n_2i = self.CppDLL.getnData2i( p_enum );   p_2i = self.CppDLL.getpData2i( p_enum )
      setattr( self, p_name, Cptr1i_to_np( p_2i, n_2i*g_nPixelFormat).reshape(n_2i,g_nPixelFormat) )

    self.CppDLL.OnDestroy.argtype = (None,)
    self.CppDLL.OnDestroy.restype = None
    self.CppDLL.OnDestroy()    

    # debug
  def Print_tabbed(self): 
    print("Va, Vb, Vtc, Vnv, Vss, Atc=")
    print( self.vm_info[0],  self.vm_info[1],# print data for  MS excel
            self.vm_info[2],  self.vm_info[3],
            self.vm_info[4],  g_mesh0_surface_area, sep=" ")  
    print("Mclad,Mcore,Mo=")
    print( self.vm_info[8],  self.vm_info[9], self.vm_info[10], sep=" ")  

