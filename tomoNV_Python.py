
# Python implementation of TomoNV algorithm
from tomoNV_io import *

"""
 Dependency:

 #pip install python==3.7.8
 #pip install numpy==1.21.4
 #pip install setuptools==60.2.0
 #pip install matplotlib==3.5.0
 #pip install pandas==1.3.5
 #pip install plotly==5.6.0
 #pip install scipy==1.7.3
 #pip install open3d==0.13.0 #use 0.10.0 in Colab
 #pip install ray==1.11.0 #not tested in Colab yet
 #pip install numba==0.55.1
 #pip install numpy-quaternion==2022.4.2 #https://quaternion.readthedocs.io/en/latest/
"""
#---------------------------------------------------------------------------
#RAY (parallelization module)  https://docs.ray.io/en/latest/ 
import ray 
from ray.util import ActorPool

def InitPARAMS(bUseRAY=True, bUseNUMBA = True,  bUseSlotPairing=True, bUseExplicitSS = True):
  global g_bUseRAY, g_nCPU, g_bUseSlotPairing, g_bUseExplicitSS, Cdll
  g_bUseRAY  = bUseRAY
  g_bUseSlotPairing = bUseSlotPairing
  g_bUseExplicitSS  = bUseExplicitSS
  import os
  os.environ["NUMBA_DISABLE_JIT"] = '0' if(bUseNUMBA) else '1'
  if(g_bUseRAY):
    import os
    g_nCPU = os.cpu_count()-1
    print('Number of CPU used by RAY:', g_nCPU)
    ray.init(num_cpus = g_nCPU)  # ini Ray parallelizatoin. takes several seconds.

@ray.remote(num_returns=1)
def getMeshInfo_RAY(param, vtx, ele, bVerbose=False):## threaded version
  return getMeshInfo(param, vtx, ele, bVerbose)

def getMeshInfo(param, vtx, ele, bVerbose=False):#non-threaded veresion
  Ray_thread_id = int(param[4])#for debug
  tomo1 = tomoNV('', bVerbose)
  tomo1.ImportMeshFromNumpy( vtx, ele)
  (tomo1.yaw, tomo1.pitch, tomo1.roll, tomo1.theta_c) = [ param[p] for p in range(0,4)]
  tomo1.Pixelize(ray_level = 2)
  tomo1.Calculate()
  if(bVerbose):
    print("ThreadID=", Ray_thread_id,", param=", FStr(param), "Mtotal=", FStr(tomo1.Mtotal))
  return (tomo1.Mo, tomo1.Mss, tomo1.Mtotal, Ray_thread_id)
#--------------------------------------------------------------
#com.geom.
def Area2( A,B,C):#assume z-component zero
  return np.cross( B-A, C-A)[2] 

def triCoord( _point, _triA, _triB, _triC):#assume z-component zero
  (x,y,z) = _point;   ax, ay, az = _triA;  bx, by, bz = _triB;  cx, cy, cz = _triC;
  area0 = Area2( _triA, _triB, _triC)
  side_1 = (x - bx) * (ay - by) - (ax - bx) * (y - by)
  side_2 = (x - cx) * (by - cy) - (bx - cx) * (y - cy)
  side_3 = (x - ax) * (cy - ay) - (cx - ax) * (y - ay)
  return (side_2 / area0, side_3/area0, side_1/area0)

import numba as nb
from numba import jit, njit

@njit( nb.types.float32[:]( nb.types.float32[:],nb.types.float32[:], nb.types.float32[:], nb.types.float32[:]))
def getBaryCoord( p, a, b, c):
  v0 = b - a; v1 = c - a; v2 = p - a
  d00 = np.dot(v0,v0); d01 = np.dot(v0,v1); d11 = np.dot(v1,v1)
  d20 = np.dot(v2,v0); d21 = np.dot(v2,v1)
  denom = d00*d11 - d01*d01
  if( np.abs(denom) > g_fMARGIN):
    v = (d11*d20 - d01*d21)/denom
    w = (d00*d21 - d01*d20)/denom
    u = 1.0 - v - w
    return np.array([u,v,w]).astype(np.float32)
  return np.array([-1.,-1.,-1.]).astype(np.float32)

def selectiveMeshCopy(_mesh0, _bMaskList):
  mesh1 = copy.deepcopy(_mesh0)
  mesh1.triangles        = o3d.utility.Vector3iVector( np.asarray( _mesh0.triangles)[ _bMaskList] )
  mesh1.triangle_normals = o3d.utility.Vector3dVector( np.asarray( _mesh0.triangle_normals)[ _bMaskList] )
  return mesh1
#--------------------------------------------------------------
#TomoNV pixel operations (Python version)

def deleteNoise(pixels, zmin, zmax):
  zcomp = pixels[:,2]
  zcomp[ zcomp < zmin ] = zmin
  zcomp[ zcomp > zmax ] = zmin
  pixels[:,2] = zcomp
  return pixels

def deleteBottomPxls( pxls):# SS points on a bottom plate are not needed.
  above_bottom = np.not_equal( pxls[:,2], 0)
  return pxls[above_bottom]


def uniqueIntList( pixels): # reduce multiple points with a same (x,y,z) into one
  pixels = pixels.astype(np.int32)
  _,idx = np.unique(pixels[:,0:3],axis=0, return_index=True)
  return pixels[idx]

@njit( nb.types.int32[:,:]())
def getEmptyPxls():
  return np.empty((0,g_nPixelFormat), dtype=np.int32)

def getEmptySlots():
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  nRow = (ymax-ymin+1) ; nCol = (xmax-xmin+1)
  return[[ getEmptyPxls() for i in range(nRow)] for j in range(nCol)] 

def slotAXPY( const_a, slot_a, const_b, slot_b):
  pxls = getEmptyPxls()
  for pxl in slot_a:
    pxls = np.append( pxls, const_a * np.array([pxl]), axis=0)
  for pxl in slot_b:
    pxls = np.append( pxls, const_b * np.array([pxl]), axis=0)
  return pxls

def slotsAXPY( const_a, slots_a, const_b, slots_b):
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  new_slots = getEmptySlots( )
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      new_slots[x][y] = slotAXPY(const_a, slots_a[x][y], const_b, slots_b[x][y])
  return new_slots
 
def getHighestPixels( pixels0, AABB2D): # [ (int64) x, y, z, nZ] 
  (xmin, ymin, xmax, ymax) = AABB2D
  buffer = createZeroPixels( xmin, ymin, xmax, ymax) 
  for pixel in pixels0:
    (pX, pY, pZ, nX, nY, nZ) = pixel
    if buffer[ pX - xmin, pY - ymin, 2] < pZ:
      buffer[  pX - xmin, pY - ymin, 2] = pZ
  pxl_map = np.reshape( np.copy( buffer), (-1, g_nPixelFormat)) #change array size for rendering
  return pxl_map

# @jit( nopython=True) ??error
def pxlsToSlots(pxls):# grouping pxls w.r.t. (x,y) coordinates
  global g_nPixelFormat
  global g_AABB2D
  (x0,y0,x1,y1) = g_AABB2D
  slots = getEmptySlots()
  for pxl in pxls:
    (x,y,z,nX,nY,nZ) = pxl
    data = np.array([[ x,y,z,nX,nY,nZ]])
    if x>= x0 and x <= x1 and y>=y0 and y<=y1:
      slots[int(x-x0)][int(y-y0)] = np.append( slots[int(x-x0)][int(y-y0)], data, axis=0)
  return slots  

#def slotsToPxls( slots, AABB2D):
def slotsToPxls( slots):
  pxls = getEmptyPxls()
  if( len(slots) ==0):
    return pxls
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      pxls = np.append( pxls, slots[x][y], axis=0)
  return pxls  

def sortSlot(slot, bDescend):
  if bDescend:
    sorting = slot[:,2].argsort()
    slot = slot[sorting[::-1]]
  else:
    slot = slot[slot[:,2].argsort()]  
  return slot

def sortSlotsByZ(slots, bDescend = True):
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      if slots[x][y].size > g_nPixelFormat:  
        slots[x][y] = sortSlot(slots[x][y], bDescend)
  return slots

def filterZNear( slot0):
  slot0 = slot0[slot0[:,2].argsort()]#sort w.r.t. z value
  slot1 = getEmptyPxls()
  last_z = -100
  for pxl in slot0:
    (x,y,z,nX,nY,nZ) = pxl
    if( z-last_z > 1 ):
      slot1 = np.append(slot1, [pxl], axis=0)
      last_z = z
  return slot1

def removeZNearPxls( slots):
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      slots[x][y] = filterZNear( slots[x][y])
  return (slots)

def getHighestPxls( slots0):
  pxls = getEmptyPxls()
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      if slots0[x][y].size > 0:
        zHighestPxl = slots0[x][y][0] #slots0 is already sorted.
        pxls = np.append( pxls, np.array([zHighestPxl]), axis=0)
  return pxls

def addZ(pxls):
  new_pxl = copy.deepcopy(pxls[0])
  new_pxl[2] = pxls.sum( axis=0)[2]
  return new_pxl

def addZs(slots0):
  pxls = getEmptyPxls()
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      if slots0[x][y].size > 0:
        new_pxl = addZ(slots0[x][y])
        pxls = np.append( pxls, np.array([new_pxl]), axis=0)
  return pxls

@njit( nb.types.float32[:]( nb.types.float32[:],nb.types.float32[:], nb.types.float32[:], nb.types.float32[:]), fastmath=True)
def bary_product( x,y,z, b_crd):
  p = x*b_crd[0] + y*b_crd[1] + z * b_crd[2]
  return p
#-------------------------------------------------------------------------------
# pixelization

#@njit(fastmath=True) 
def triPixel( tri0):#tri1 = [v0, v1, v2]
  global g_fMARGIN, g_HALF_VOXEL_SIZE
  tri2d = np.copy(tri0)
  vtx = tri2d[0:3] + g_fMARGIN;nrm=tri2d[3:6]
  vtx_min = [np.min( vtx[:,0]), np.min( vtx[:,1]), np.min( vtx[:,2])]
  vtx_max = [np.max( vtx[:,0]), np.max( vtx[:,1]), np.max( vtx[:,2])]
  (x0,y0,z0) = list(map(int, vtx_min))
  (x1,y1,z1) = list(map(int, vtx_max))
  vtx[:,2]=0 #move pixels to xy plan
  x = np.arange( x0, x1+2).astype(np.float32)
  y = np.arange( y0, y1+2).astype(np.float32)
  x_tile = np.repeat(x, len(y))
  y_tile = np.repeat(y, len(x)).reshape(-1,len(x)).T.flatten()
  coord = np.dstack( (x_tile, y_tile)).reshape( len(x), len(y), 2)
  pxls  = getEmptyPxls()
  for i in range(len(x)):
    for j in range(len(y)):
      v_center =  np.array( [coord[i,j,0] + g_fMARGIN + g_HALF_VOXEL_SIZE,coord[i,j,1]+ g_fMARGIN + g_HALF_VOXEL_SIZE, g_HALF_VOXEL_SIZE] ).astype(np.float32)
      b_crd = getBaryCoord(v_center, vtx[0], vtx[1], vtx[2])
      if(b_crd[0] >= -g_fMARGIN and b_crd[1] >= -g_fMARGIN and b_crd[2] <= 1.+ g_fMARGIN and b_crd[0] + b_crd[1] <= 1.+g_fMARGIN):
        xyz = bary_product( tri0[0], tri0[1], tri0[2], b_crd) #pxl coordinates 
        N   = bary_product( nrm[0],  nrm[1],  nrm[2], b_crd) #pxl coordinates 
        new_pxl = np.array([[ 
          int(round(xyz[0],1)), 
          int(round(xyz[1],1)), 
          int(round(xyz[2],1)),
          int(round(N[0]*g_fNORMALFACTOR,1)),
          int(round(N[1]*g_fNORMALFACTOR,1)),
          int(round(N[2]*g_fNORMALFACTOR,1))
          ]]) 
        pxls = np.append( pxls, new_pxl, axis=0)
  if pxls.size == 0:#very small triangle.
    xyz = (np.asarray(vtx_min) + np.asarray(vtx_max)) * 0.5
    N = copy.deepcopy(nrm[0]*g_fNORMALFACTOR)
    new_pxl = np.array([[    xyz[0], xyz[1], xyz[2], N[0],N[1], N[2]    ]]).astype(np.int32)
    pxls = np.append( pxls, new_pxl, axis=0)    
  return pxls #  [(integers) x, y, z, nx, ny, nz]

@ray.remote(num_returns=1)
def _pixelizeMesh_RAY( vtx, ele, vtxNrm): #parellellized version
  return _pixelizeMesh( vtx, ele, vtxNrm) 

def _pixelizeMesh( vtx, ele, vtxNrm): #default version
  pxls = getEmptyPxls()
  for e in ele:
    tri = np.array( [ 
        vtx[e[0]],    vtx[e[1]],     vtx[e[2]],
        vtxNrm[e[0]], vtxNrm[e[1]],  vtxNrm[e[2]]     ]) # [v0, v1, v2]
    new_pxls = triPixel( tri)
    if new_pxls.size > 0:
      pxls = np.append( pxls, new_pxls, axis=0)
  return pxls


def pixelizeMesh( _mesh0, ray_level): # non-parellel version.
  global g_bUseRAY
  vtx = np.asarray( _mesh0.vertices).astype(np.float32)
  ele = np.asarray( _mesh0.triangles).astype(np.int32)
  nrm = np.asarray( _mesh0.vertex_normals).astype(np.float32)#for NV pixels. use per-pixel Gouraud normal, to erase overlapping pixels exactly.
  if(g_bUseRAY and ray_level==1):
    vtx_id = ray.put( vtx); ele_id = ray.put( ele); nrm_id = ray.put( nrm) 
    return ray.get(      _pixelizeMesh_RAY.remote( vtx_id, ele_id, nrm_id))
  else:
    return _pixelizeMesh( vtx, ele, nrm)

#-------------------------------------------------------------------------------
# slot pairing

def getShadowCastor( be_pxls, theta_c):#explicit + implicit
  nrm    = be_pxls[:,3:6].astype(np.float32) / g_fNORMALFACTOR
  bSS    = np.dot( nrm, zaxis) < - np.sin( theta_c )
  bNV    = np.logical_not( bSS)
  ssb_pixels = be_pxls[bSS]
  nvb_pixels = be_pxls[bNV]
  return (ssb_pixels, nvb_pixels)

def getNVA( nvb, al_slot):#return alpha point next to nvb in Z-direction
  sorting = al_slot[:,2].argsort()
  al_slot = al_slot[sorting[::-1]]#inverse sorting w.r.t. z value. https://codetorial.net/tips_and_examples/numpy_argsort.html
  z0 = nvb[2]
  for al_pxl in al_slot:
    if( al_pxl[2] <= z0):#//등호 꼭 필요함.
      return al_pxl
  bottom = nvb; bottom[2] = 0
  return bottom # return bottom plate

# @njit( nb.types.UniTuple(nb.types.int32[:,:], 2)( nb.types.int32[:,:], nb.types.int32[:,:]) )
def _matchPairNumber( al_slot, be_slot):
  if al_slot.size < be_slot.size:
    sorting = be_slot[:,2].argsort()
    be_slot = be_slot[sorting]
    be_slot = be_slot[1:al_slot.size:1]
  elif al_slot.size > be_slot.size:
    sorting = al_slot[:,2].argsort()
    al_slot = al_slot[sorting[::-1]]
    al_slot = al_slot[1:be_slot.size:1]
  return (al_slot, be_slot)

@njit( nb.types.int32[:,:]( nb.types.int32[:,:], nb.types.int32[:,:]) )
def _matchPairWRT( pxls, ref_pxls):
  if pxls.size > ref_pxls.size:
    pxls = pxls[1:ref_pxls.size:1]
  return pxls

@njit( nb.types.boolean( nb.types.int32[:,:], nb.types.int32, nb.types.int32) )
def _hasPxlBetween(slot, z_low, z_high):#assume SLOT is already sorted w.r.t. z-coordinate 
  for pxl in slot:
    if (pxl[2] <= z_high and pxl[2] >= z_low) or (pxl[2] <= z_low and pxl[2] >= z_high):
      return True 
  return False

# @njit( nb.types.int32[:,:]( nb.types.int32[:,:], nb.types.int32[:,:]) ) #??error
def _matchAlternationWRT( slot, ref_slot):
  (nPxl,_) = np.shape(slot)
  if nPxl <= 1:
    return slot
  new_slot = getEmptyPxls()
  for s in range(nPxl-1):
    pxl0 = slot[s]; pxl1 = slot[s+1]
    if _hasPxlBetween( ref_slot, pxl0[2],pxl1[2]):
      new_slot = np.append( new_slot, np.array([pxl0]), axis=0)
  new_pxl = np.array([ slot[nPxl-1] ])
  new_slot = np.append( new_slot, new_pxl, axis=0)# dd last(highest) point to preserve volume
  return new_slot

def matchAlBeSlots( al_slots, be_slots):
  #pair matcing for noise suppression (Aassume that number of alpha and beta pixels are same for closed volume object)
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      al_slots[x][y] = _matchAlternationWRT( al_slots[x][y], be_slots[x][y])
      be_slots[x][y] = _matchAlternationWRT( be_slots[x][y], al_slots[x][y])
      (al_slots[x][y], be_slots[x][y]) = _matchPairNumber( al_slots[x][y], be_slots[x][y])
  return (al_slots, be_slots)

def matchPairsWRT( targt_pxls, ref_pxls, AABB2D):
  ref_slots = pxlsToSlots( ref_pxls)
  targt_slots = pxlsToSlots( targt_pxls)
  new_pxls = getEmptyPxls()
  (xmin, ymin, xmax, ymax) = AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      reduced_nvs = _matchPairWRT( targt_slots[x][y], ref_slots[x][y])
      new_pxls = np.append( new_pxls, reduced_nvs, axis=0)
  return ( new_pxls)

def matchSlotsWRT( slots, ref_slots):
  global g_AABB2D
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      slots[x][y] = _matchPairWRT( slots[x][y], ref_slots[x][y])
  return ( slots)

def createTCSlots( al_slots):
  global g_AABB2D
  TC_slots = getEmptySlots()
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      if al_slots[x][y].size > 0:
        TC_slots[x][y] = np.append( TC_slots[x][y], np.array([al_slots[x][y][0]]), axis=0)
  return TC_slots

def createShadowCastor( be_slots, theta_c):#find explicit/implicit pxls at a time
  global g_AABB2D
  ssb_slots = getEmptySlots()
  nvb_slots = getEmptySlots()
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      if be_slots[x][y].size > 0:
        (ssb_slots[x][y], nvb_slots[x][y]) = getShadowCastor( be_slots[x][y], theta_c)
  return (ssb_slots, nvb_slots)

def getShadowAcceptor( ssb_pxl, ssa_pxl):#ssb_pxl, ssa_pxl have same (x,y) slot. only differs in Z coordinate
  ss_pxls   = getEmptyPxls()
  for z_offset in range(ssb_pxl[2] - ssa_pxl[2]+1):
    new_ss  = copy.deepcopy(ssb_pxl)
    new_ss[2] -= z_offset
    ss_pxls = np.append( ss_pxls, np.array([new_ss]), axis=0)
  return ss_pxls

def createShadowAcceptorSlots( NVB_slots, al_slots):
  global g_AABB2D
  NVA_slots = getEmptySlots()
  ss_pxls   = getEmptyPxls()#for rendering of explicit SS
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):
      for nvb in NVB_slots[x][y]:
        if nvb[2] > 0:
          new_nbA = getNVA( copy.deepcopy(nvb),  al_slots[x-xmin][y-ymin])
          NVA_slots[x][y] = np.append( NVA_slots[x][y], np.array([new_nbA]), axis=0)  
          new_ss_pxls  = getShadowAcceptor( nvb, new_nbA)
          if(new_ss_pxls.size > 0):
            ss_pxls = np.append(ss_pxls,np.array(new_ss_pxls), axis=0)  
  return (NVA_slots, ss_pxls)

def createSSSlots_Implicit( al_slots, be_slots, tc_slots, nvb_slots, nva_slots):
  global g_AABB2D
  ss_slots = getEmptySlots()
  (xmin, ymin, xmax, ymax) = g_AABB2D
  for x in range(xmax-xmin+1):
    for y in range(ymax-ymin+1):

      if( x== 7 and y == 3):
        debug = 1

      alpha_minus = np.copy(al_slots[x][y]);  alpha_minus[:,2] *= -1
      nvB_minus   = np.copy(nvb_slots[x][y]); nvB_minus[:,2] *= -1
      tc          = np.copy(tc_slots[x][y])
      alpha_minus[:,g_nPixelFormat-1] = 0;  nvB_minus[:,g_nPixelFormat-1] = 0;  tc[:,g_nPixelFormat-1] = 0  # we need only Nz of beta pixels
      ss = alpha_minus
      ss = np.append( ss, be_slots[x][y],   axis=0)
      ss = np.append( ss, tc,               axis=0)
      ss = np.append( ss, nvB_minus,        axis=0)
      ss = np.append( ss, nva_slots[x][y],  axis=0)
      ss_slots[x][y] = ss
  return ss_slots

#--------------------------------------------------------------
# searaching optimal orientation.

def SearchYP( tomo0, yaw_range, pitch_range, bVerbose=False, bShowGraph=False):
  vtx = np.array(tomo0.mesh0.vertices)
  ele = np.array(tomo0.mesh0.triangles)
  (params, sizeY, sizeP, sizeYP) = paramsYP( yaw_range, pitch_range, tomo0.theta_c)

  if(bVerbose):
    print('(yaw,pitch,roll) list to search =\n', FStr( toDegree(params[0:3])) )

  #multi-threading for nBatch groups..
  tomo0.Mo3D     = np.zeros( sizeYP, dtype=np.float32) 
  tomo0.Mss3D    = np.zeros( sizeYP, dtype=np.float32) 
  tomo0.Mtotal3D = np.zeros( sizeYP, dtype=np.float32) 
  global g_bUseRAY
  if(g_bUseRAY):
    vtx_id = ray.put( vtx);    ele_id = ray.put( ele)
    ray_IDs = [ getMeshInfo_RAY.remote( p, vtx_id, ele_id, bVerbose )  for p in params ]
    while len(ray_IDs) > 0: # https://docs.ray.io/en/latest/ray-core/examples/tips-for-first-time.html#tip-4-pipeline-data-processing
      finished_IDs, ray_IDs = ray.wait(ray_IDs)
      result = ray.get(finished_IDs[0])
      p = result[3]
      (tomo0.Mo3D[p], tomo0.Mss3D[p],  tomo0.Mtotal3D[p], thread_id) = result
  else:
    for p in range(sizeYP):
      (tomo0.Mo3D[p], tomo0.Mss3D[p], tomo0.Mtotal3D[p], thread_id) = getMeshInfo( params[p], vtx, ele, bVerbose)

  #print optimal orientation info.
  global g_nOptimalsToDisplay
  (tomo0.Optimals,tomo0.Worsts) = findOptimals(params, tomo0.Mss3D, g_nOptimalsToDisplay) #Caution: data should be 1D array.
  tomo0.Mo3D     = np.reshape( tomo0.Mo3D,     (sizeY,sizeP))
  tomo0.Mss3D    = np.reshape( tomo0.Mss3D,    (sizeY,sizeP))
  tomo0.Mtotal3D = np.reshape( tomo0.Mtotal3D, (sizeY,sizeP))

  #show search result graphically
  if(bShowGraph):
    Plot3D( tomo0.mesh0, yaw_range, pitch_range, tomo0.Mss3D, tomo0.Optimals, tomo0.Worsts)

  #if(bVerbose):
  print("Mo3D=",      FStr( tomo0.Mo3D))
  print("Mss3D=",     FStr( tomo0.Mss3D))
  print("Mtotal3D=",  FStr( tomo0.Mtotal3D))
  np.savetxt("Mo3D.txt", tomo0.Mo3D, delimiter='\n', fmt='%.2f')
  np.savetxt("Mss3D.txt", tomo0.Mss3D, delimiter='\n', fmt='%.2f')
  np.savetxt("Mtotal3D.txt", tomo0.Mtotal3D, delimiter='\n', fmt='%.2f')
  
  return( tomo0.Mo3D, tomo0.Mss3D, tomo0.Mtotal3D)
#end of def SearchYP()

def paramsYPR( yaw_range, pitch_range, roll_range, theta_c):
  params = np.empty( [0,5], dtype=np.float32) 
  sizeY = yaw_range.size ; sizeP = pitch_range.size ; sizeR = roll_range.size; sizeYPR = sizeY * sizeP * sizeR
  thread_id = 0#RAY deug
  for yaw in yaw_range:
    for pitch in pitch_range:
      for roll in  roll_range:
        params = np.append( params, np.array([[yaw,pitch,roll,theta_c, thread_id+ g_fMARGIN]]), axis=0)
        thread_id = thread_id + 1  
  return ( params, sizeY, sizeP, sizeR, sizeYPR)


def SearchYPR( tomo0, yaw_range, pitch_range, roll_range, bVerbose=False, bShowGraph=False):
  (params, sizeY, sizeP, sizeR, sizeYPR) = paramsYPR(yaw_range, pitch_range, roll_range, tomo0.theta_c)
  if(bVerbose):
    print('(yaw,pitch,roll,theta_c=\n', FStr( toDegree(params)) )

  #iteration
  vtx = np.array(tomo0.mesh0.vertices)
  ele = np.array(tomo0.mesh0.triangles)
  tomo0.Mo4D     = np.zeros( sizeYPR, dtype=np.float32) 
  tomo0.Mss4D    = np.zeros( sizeYPR, dtype=np.float32) 
  tomo0.Mtotal4D = np.zeros( sizeYPR, dtype=np.float32) 
  global g_bUseRAY
  if(g_bUseRAY):#multithreading using RAY component
    vtx_id = ray.put( vtx)
    ele_id = ray.put( ele)
    ray_results = ray.get( [ getMeshInfo_RAY.remote( p, vtx_id, ele_id, bVerbose=False )  for p in params] )
    for p in range(sizeYPR):
      tomo0.Mo4D[p]      = ray_results[p][0]
      tomo0.Mss4D[p]     = ray_results[p][1]
      tomo0.Mtotal4D[p]  = ray_results[p][2]
    if(bVerbose):
      print( "Ray(multithread)_results=" , FStr(ray_results))
  else:#non-threaded version
    for p in range(sizeYPR):
       (tomo0.Mo4D[p], tomo0.Mss4D[p], tomo0.Mtotal4D[p]) = getMeshInfo( params[p], vtx, ele, bVerbose=False)

  #print optimal orientation info.
  print("Mo4D=",      FStr( tomo0.Mo4D))
  print("Mss4D=",     FStr( tomo0.Mss4D))
  print("Mtotal4D=",  FStr( tomo0.Mtotal4D))
  np.savetxt("Mo4D.txt", tomo0.Mo4D, delimiter='\n', fmt='%.2f')
  np.savetxt("Mss4D.txt", tomo0.Mss4D, delimiter='\n', fmt='%.2f')
  np.savetxt("Mtotal4D.txt", tomo0.Mtotal4D, delimiter='\n', fmt='%.2f')

  (tomo0.Optimals,tomo0.Worsts) = findOptimals(params, tomo0.Mtotal4D) #Caution: data should be 1D array.
  tomo0.Mo4D     = np.reshape( tomo0.Mo4D,     (sizeY,sizeP,sizeR))
  tomo0.Mss4D    = np.reshape( tomo0.Mss4D,    (sizeY,sizeP,sizeR))
  tomo0.Mtotal4D = np.reshape( tomo0.Mtotal4D, (sizeY,sizeP,sizeR))

  #draw graph
  if(bShowGraph):
    (nOpt, aaa) = tomo0.Optimals.shape
    gs_kw = dict(width_ratios=[2., 1,1.,1.,1.,1.], height_ratios=[1])
    (fig, axes) = plt.subplots(nrows=1, ncols=nOpt+1, figsize=(4*nOpt,4),
      gridspec_kw=gs_kw,constrained_layout=True, subplot_kw={"projection":"3d"}) #  #https://towardsdatascience.com/creating-custom-plotting-functions-with-matplotlib-1f4b8eba6aa1
    (X, Y, Z) = np.meshgrid( yaw_range, pitch_range, roll_range)
    graphYPR_Plotly(  X, Y, Z, tomo0.Mtotal4D, tomo0.Optimals)
    graphYPR( X, Y, Z, tomo0.Mtotal4D, tomo0.Optimals, ax = axes[0])
    [ drawOptimals( tomo0, opt, 'optimal', ax = axes[opt+1]) for opt in range(nOpt) ]
    plt.show()
    fig.savefig( "Mtotal4D.png", dpi=200)

  return( tomo0.Mo4D, tomo0.Mss4D, tomo0.Mtotal4D)
  #end of def SearchYPR()

#==================================================================================================
# Python version TomoNV main class definition

class tomoNV:
  #default values
  (yaw, pitch, roll, theta_c)   = (0. , 0., 0., toRadian(60.)) #in radian
  wall_thickness = 0.8 # [mm]
  PLA_density    = 0.00121 # density of PLA filament, [g/mm^3]
  Fclad = 1.0 # fill ratio of cladding, always 1.0
  Fcore = 0.15 # fill ratio of core, (0~1.0)
  Fss   = 0.2 # fill ratio of support structure, (0~1.0)
  Css   = 1.2**3 # correction constant for Mss. obsolete. replaced by getWeightZSum[]
  Fclad = 1.
  bVerbose = True #print intermediate messages  
  ray_level = 1 #1 = multithreading only in pixelizeMesh(), 2= multithreading of overall steps(For YP search)
  bUseSlotPairing = True
  FileName = ""
  al_pxls=[]; be_pxls=[];  NVB_pxls=[]; NVA_pxls=[]; Vss_pxls=[]; TC_pxls=[]; Vo_pxls=[]; SS_pxls=[]; SSB_pxls=[]; SSA_pxls=[]
  AABB2D = (0,0,0,0)

  def __init__(self, filename, bVerbose=False):
    global g_bUseSlotPairing, g_bUseExplicitSS
    self.bUseExplicitSS  = g_bUseExplicitSS
    self.bUseSlotPairing = g_bUseSlotPairing
    self.bVerbose = bVerbose
    for pxls in g_PixelVarNames:
      setattr(self, pxls, getEmptyPxls())
    #Load file
    import os.path    
    self.mesh0 = o3d.geometry.TriangleMesh()
    if( os.path.isfile(filename) ):
      self.FileName = filename
      self.mesh0 = o3d.io.read_triangle_mesh(filename)
      self.mesh0 = self.mesh0.translate( self.mesh0.get_center() * -1.)
      self.mesh0_surface_area = o3d.geometry.TriangleMesh.get_surface_area( self.mesh0)
      # (self.mesh0_min, self.mesh0_max,_) = getBoundary( self.mesh0.vertices)

      if(self.bVerbose): 
        print('Mesh data loaded:' + filename)

  def ImportMeshFromNumpy( self, vtx, ele):#for Get4D_Parallel[]
    self.mesh0.vertices  = o3d.utility.Vector3dVector( copy.deepcopy(vtx))
    self.mesh0.triangles = o3d.utility.Vector3iVector( copy.deepcopy(ele))
    self.mesh0_surface_area = o3d.geometry.TriangleMesh.get_surface_area( self.mesh0)
    (self.mesh0_min, self.mesh0_max,_) = getBoundary( self.mesh0.vertices)

  def Pixelize_Step1(self):# Step 1. rotate and move onto bottom plate
    self.mesh1 = copy.deepcopy( self.mesh0)

    qn    = Rotation.from_euler('xyz', [[self.yaw, self.pitch, self.roll]], degrees=False)
    vtx1      = np.asarray( self.mesh1.vertices)
    vtx1     = qn.apply(vtx1)
    vtx1_min = np.min( vtx1, axis=0)
    vtx1  -= vtx1_min
    self.mesh1.vertices = o3d.utility.Vector3dVector( vtx1)

    # self.mesh1.compute_triangle_normals() # Caution: compute_triangle_normals() gives wrong normal.
    self.mesh1.compute_vertex_normals()
    (self.mesh1_min, self.mesh1_max,_) = getBoundary( self.mesh1.vertices)
    (self.x1,self.y1,self.z1) = list(map(int, self.mesh1_max)) # slot size. integer.
    (self.x0,self.y0,self.z0) = list(map(int, self.mesh1_min)) # slot size. integer.
    
  def Pixelize_Step2(self):# Step 2. find alpha, beta triangles
    self.tri_nrm = np.asarray( self.mesh1.triangle_normals)
    self.bAlpha = np.dot( self.tri_nrm, zaxis) >  g_fMARGIN * 10.
    self.bBeta  = np.dot( self.tri_nrm, zaxis) < -g_fMARGIN * 10.
    self.al_mesh = selectiveMeshCopy( self.mesh1, self.bAlpha)
    self.be_mesh = selectiveMeshCopy( self.mesh1, self.bBeta)

  def Pixelize_Step3(self):# Step 3. make triangle groups. time consuming.
    self.al_pxls  = pixelizeMesh( self.al_mesh, self.ray_level)
    self.be_pxls  = pixelizeMesh( self.be_mesh, self.ray_level)

  def Pixelize_Step4(self):#step 4. make pixel groups
    global g_AABB2D
    g_AABB2D = (self.x0, self.y0, self.x1, self.y1)
    al_slots  = pxlsToSlots( self.al_pxls)
    be_slots  = pxlsToSlots( self.be_pxls)

    #Pair-matching to reduce noise pixls
    if self.bUseSlotPairing:
      al_slots = sortSlotsByZ(removeZNearPxls( al_slots))
      be_slots = sortSlotsByZ(removeZNearPxls( be_slots))
      (al_slots, be_slots) = matchAlBeSlots( al_slots, be_slots)#main step      

    Vo_slots  = slotsAXPY( 1., al_slots, -1., be_slots)
    TC_slots  = createTCSlots( al_slots)#for implicit SS

    if(self.bUseExplicitSS):
      #explicit form of s.s. pixels
      (SSB_slots,_) = createShadowCastor( be_slots, self.theta_c)
      SSB_slots     = removeZNearPxls( SSB_slots)
      SSB_slots     = matchSlotsWRT( SSB_slots, be_slots)# No. of SSB cannot be greater than no. of beta.
      (SSA_slots, self.SS_pxls)  = createShadowAcceptorSlots( copy.deepcopy(SSB_slots), al_slots)
      Vss_slots  = slotsAXPY( 1., SSB_slots, -1., SSA_slots)
    else:
      #implicit form of s.s. pixels
      (_,NVB_slots) = createShadowCastor( be_slots, self.theta_c)
      NVB_slots     = removeZNearPxls( NVB_slots)
      NVB_slots     = matchSlotsWRT( NVB_slots, be_slots)# No. of NVB cannot be greater than no. of beta.
      (NVA_slots, _)  = createShadowAcceptorSlots( copy.deepcopy(NVB_slots), al_slots)
      Vss_slots  = createSSSlots_Implicit( al_slots, be_slots, TC_slots, NVB_slots, NVA_slots)

    #slots to pxls
    self.al_pxls  = slotsToPxls( al_slots)
    self.be_pxls  = slotsToPxls( be_slots)
    self.Vo_pxls  = addZs( Vo_slots)
    self.TC_pxls  = uniqueIntList( slotsToPxls( TC_slots))
    self.Vss_pxls = addZs( Vss_slots) 

    #these are for debugging (not necessary)
    if(self.bUseExplicitSS):
      self.SSB_pxls = slotsToPxls(SSB_slots)
      self.SSA_pxls = slotsToPxls(SSA_slots)
    else:
      self.NVB_pxls = slotsToPxls(NVB_slots)
      self.NVA_pxls = slotsToPxls(NVA_slots)

    self.Vss_pxls  = deleteBottomPxls( self.Vss_pxls) #for rendering
    self.Vo_pxls   = deleteBottomPxls( self.Vo_pxls) #for rendering
    
    if not self.bUseSlotPairing:
      self.Vo_pxls  = deleteNoise( self.Vo_pxls,  self.z0, self.z1)
      self.Vss_pxls = deleteNoise( self.Vss_pxls, self.z0, self.z1)
  
  def Pixelize(self, ray_level=1):#partially parallelize version 
    self.ray_level = ray_level
    self.Pixelize_Step1()
    self.Pixelize_Step2()
    self.Pixelize_Step3()#Use parallel computation in pixelizeMeshRAY[] function
    self.Pixelize_Step4()
  
    if(self.bVerbose): 
      self.Calculate() #to print Mtotal value
      print('Finding pixels at [', 
          FStr(toDegree(self.yaw)), 
          FStr(toDegree(self.pitch)),
          FStr(toDegree(self.roll)),']\n',
          '          -> Mo=', FStr(self.Mo),     "[g],",
                      'Mss=', FStr(self.Mss),    "[g]",
                   'Mtotal=', FStr(self.Mtotal), "[g]")
    
  def Calculate(self):
    # calculate volume
    self.Va  = self.al_pxls.sum( axis=0)[2] * (g_voxel_size*g_voxel_size)
    self.Vb  = self.be_pxls.sum( axis=0)[2] * (g_voxel_size*g_voxel_size)
    self.Vtc = self.TC_pxls.sum( axis=0)[2] * (g_voxel_size*g_voxel_size)
    self.Vnv = (self.NVB_pxls.sum( axis=0)[2] - self.NVA_pxls.sum( axis=0)[2]) * (g_voxel_size*g_voxel_size)
    self.Vo  = self.Vo_pxls.sum( axis=0)[2] * (g_voxel_size*g_voxel_size)#=  Va - Vb

    self.Vss = self.Vss_pxls.sum( axis=0)[2] * (g_voxel_size*g_voxel_size) #self.Vss = -self.Va + self.Vb + self.Vtc - self.Vnv
 
    self.Vclad = (self.mesh0_surface_area) * self.wall_thickness 
    self.Vcore = self.Vo - self.Vclad
    
    # calculate mass
    self.Mclad = self.Vclad * self.Fclad * self.PLA_density
    self.Mcore = self.Vcore * self.Fcore * self.PLA_density
    self.Mss = self.Vss * self.Fss * self.PLA_density * self.Css
    self.Mo = self.Mcore + self.Mclad
    self.Mtotal = self.Mo + self.Mss
     
  def Print(self):
    #print out
    print("-------Volume info.----")
    print('Va=', self.Va, ', Vb=',  self.Vb, ', Vtc=',   self.Vtc,  ', Vnv=',   self.Vnv)
    print('Vo=', self.Vo, ', Vss=', self.Vss,', Vclad=', self.Vclad,', Vcore=', self.Vcore)
    print("-------Mass info.------")
    print("Mcore=", self.Mcore,", Mclad=", self.Mclad)
    print("Mo=", self.Mo, ", Mss=", self.Mss, ", Mtotal=", self.Mtotal)

  def Print_tabbed(self, X=-1,Y=-1): 
    # print("[",self.Filename,"]")
    print( self.Va,  self.Vb,# print data for  MS excel
           self.Vtc, self.Vnv,
           self.Vss, self.mesh0_surface_area, sep=" ")  
           
  def FigureTitle(self):
    return str(FStr(toDegree([self.yaw, self.pitch, self.roll]),precision=0))

# end of main class definition ===============================================================
