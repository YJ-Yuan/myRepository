# Magnetic field unit in this code is mT
# Assume each loop coil is square for calculation simplicity
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def crossProduct(a,b):
	x = a[1]*b[2]-a[2]*b[1]
	y = a[2]*b[0]-a[0]*b[2]
	z = a[0]*b[1]-a[1]*b[0]
	return (x,y,z)

#loop(xOrigin,Radius,directionSign,meshNumber)
#loop().MField(positionForField,Current)
class loop:
	def __init__(self,b,e,mesh):
		begin = [0,0,0]
		end = [0,0,0]
		begin[0] = b[0]-0.0427/2
		end[0] = e[0]-0.05
		begin[1] = b[1] -0.0427/2
		begin[2] = b[2] -0.0427/2
		end[0] = e[0]-0.0427/2
		end[1] = e[1] -0.0427/2
		end[2] = e[2] -0.0427/2
		self.points = []
		dx = (end[0]-begin[0])/mesh
		dy = (end[1]-begin[1])/mesh
		dz = (end[2]-begin[2])/mesh
		self.dl = sqrt(dx**2+dy**2+dz**2)
		for i in range(mesh):
			x = begin[0]+i*dx
			y = begin[1]+i*dy
			z = begin[2]+i*dz
			position = (x,y,z)
			direction = (dx/self.dl,dy/self.dl,dz/self.dl)
			point = [position,direction]
			self.points.append(point)


	# dB = 1/(10*d^2) * Idl x r
	def MField(self,pos,I):
		intensityX = 0
		intensityY = 0
		intensityZ = 0
		for point in self.points:
			position = point[0] #position = (x,y,z)
			direction = point[1] #direction = (x,y,z)
			dlVec = tuple(self.dl*i for i in direction) #dl for cross product
			r = tuple(pos[i]-position[i] for i in range(3))

			absr = 0
			for i in range(3):
				absr += r[i]**2
			absr = sqrt(absr)

			tmp = crossProduct(dlVec, r)
			intensityX += I*tmp[0]/1000/absr**3
			intensityY += I*tmp[1]/1000/absr**3
			intensityZ += I*tmp[2]/1000/absr**3
		pass
		return (intensityX,intensityY,intensityZ)
	pass
pass

class solenoid:
	# centerX is the center of X, R is the smallest radius of all solenoids
	def __init__(self,mesh,I):
		self.I = I
		self.loops = []
		self.loops.append(loop((0,0,0),(0,0,0.0427),mesh))
		self.loops.append(loop((0,0,0),(0,0.0427,0),mesh))
		self.loops.append(loop((0,0.0427,0.0427),(0,0.0427,0),mesh))
		self.loops.append(loop((0,0.0427,0.0427),(0,0,0.0427),mesh))
		self.loops.append(loop((0.0427,0,0.0427),(0.0427,0.0427,0.0427),mesh))
		self.loops.append(loop((0.0427,0.0427,0),(0.0427,0.0427,0.0427),mesh))
		self.loops.append(loop((0.0427,0.0427,0),(0.0427,0,0),mesh))
		self.loops.append(loop((0.0427,0,0.0427),(0.0427,0,0),mesh))

		self.loops.append(loop((0.0427,0.0427,0.0427),(0,0.0427,0.0427),mesh))
		self.loops.append(loop((0,0,0.0427),(0.0427,0,0.0427),mesh))
		self.loops.append(loop((0.0427,0,0),(0,0,0),mesh))
		self.loops.append(loop((0,0.0427,0),(0.0427,0.0427,0),mesh))
		self.loops.append(loop((0.0427,0.0427,0.0427),(0,0.0427,0.0427),mesh))
		self.loops.append(loop((0,0,0.0427),(0.0427,0,0.0427),mesh))
		self.loops.append(loop((0.0427,0,0),(0,0,0),mesh))
		self.loops.append(loop((0,0.0427,0),(0.0427,0.0427,0),mesh))

		

	def SumField(self,pos):
		intensity = [0,0,0]
		for l in self.loops:
			tmp = l.MField(pos,self.I)
			intensity[0] += tmp[0]
			intensity[1] += tmp[1]
			intensity[2] += tmp[2]
		return intensity

	def meshXZ(self,xMax,zMax):
		x0 = linspace(-xMax/100,xMax/100,30)
		z0 = linspace(-zMax/100,zMax/100,30)
		x, z = meshgrid(x0,z0,indexing='ij')
		[intensityX,intensityY,intensityZ] = self.SumField((x,0,z))
		# fig,ax = plt.subplots()
		# ax.quiver(x,z,intensityX,intensityZ)
		inten = sqrt(intensityX**2+intensityY**2+intensityZ**2)
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.set_xlabel('X(cm)')
		ax.set_ylabel('Z(cm)')
		ax.set_zlabel('Intensity(G)')
		ax.set_title('B on XZ plane')

		# Plot the surface.
		surf = ax.plot_surface(x*100, z*100, inten, cmap=cm.coolwarm,
		                       linewidth=0, antialiased=False)
		fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()

	def meshYZ(self,yMax,zMax,x=0):
		y0 = linspace(-yMax/100,yMax/100,30)
		z0 = linspace(-zMax/100,zMax/100,30)
		y, z = meshgrid(y0,z0,indexing='ij')
		[intensityX,intensityY,intensityZ] = self.SumField((x,y,z))
		# fig,ax = plt.subplots()
		# ax.quiver(x,z,intensityX,intensityZ)
		inten = sqrt(intensityX**2+intensityY**2+intensityZ**2)
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.set_xlabel('Y(cm)')
		ax.set_ylabel('Z(cm)')
		ax.set_zlabel('Intensity(G)')
		ax.set_title('B on YZ plane')

		# Plot the surface.
		surf = ax.plot_surface(y*100, z*100, inten, cmap=cm.coolwarm,
		                       linewidth=0, antialiased=False)
		fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()
		pass

	def gradientX(self,xMax):
		x = linspace(-xMax/100,xMax/100,300)
		intensity = self.SumField((x,0,0))
		for i in range(3):
			intensity[i] = gradient(intensity[i],2)
		plt.plot(x*100,intensity[0],label="Bx")
		plt.plot(x*100,intensity[1],label="By")
		plt.plot(x*100,intensity[2],label="Bz")
		plt.xlabel('X(cm)')
		plt.ylabel('∂B/∂X(G/cm)')
		plt.title('Gradient of B vs position of X')
		plt.legend(loc=2)
		plt.show()
	def gradientZ(self,zMax):
		z = linspace(-zMax/100,zMax/100,300)
		intensity = self.SumField((0,0,z))
		for i in range(3):
			intensity[i] = gradient(intensity[i],2)
		plt.plot(z*100,intensity[0],label="Bx")
		plt.plot(z*100,intensity[1],label="By")
		plt.plot(z*100,intensity[2],label="Bz")
		plt.xlabel('Z(cm)')
		plt.ylabel('∂B/∂X(G/cm)')
		plt.title('Gradient of B vs position of Z')
		plt.legend(loc=2)
		plt.show()
	def BX(self,xMax):
		x = linspace(-xMax/100,xMax/100,300)
		intensity = self.SumField((x,0,0))
		plt.plot(x*100,intensity[0],label="Bx")
		plt.plot(x*100,intensity[1],label="By")
		plt.plot(x*100,intensity[2],label="Bz")
		plt.xlabel('X(cm)')
		plt.ylabel('B(G)')
		plt.title('B vs position of X')
		plt.legend(loc=2)
		plt.show()
	def BZ(self,zMax):
		z = linspace(-zMax/100,zMax/100,300)
		intensity = self.SumField((0,0,z))
		plt.plot(z*100,intensity[0],label="Bx")
		plt.plot(z*100,intensity[1],label="By")
		plt.plot(z*100,intensity[2],label="Bz")
		plt.xlabel('Z(cm)')
		plt.ylabel('B(G)')
		plt.title('B vs position of Z')
		plt.legend(loc=2)
		plt.show()

#centerX,R,thickness,xNum,rNum,mesh,sign,I
l = solenoid(2000,70)
l.meshXZ(10, 10)
# l.meshXZ(1.5, 1.5)
# l.meshYZ(1.5, 1.5)
# l.gradientX(3)
# l.gradientZ(2)
# l.BX(3)
# l.BZ(2)

