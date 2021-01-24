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
	def __init__(self,x,R,sign,mesh):
		self.points = []
		dθ = 2*pi / mesh
		self.dl = R * dθ
		for i in range(mesh): #range(x) is 0,1,2...x-1
			θ = dθ * i
			position = (x,R*cos(θ),R*sin(θ))
			if sign == 1:
				direction = (0,-sin(θ),cos(θ))
			elif sign == -1:
				direction = (0,sin(θ),-cos(θ))
			else:
				print("Direction error")
			point = [position,direction]
			self.points.append(point)
		pass
	pass

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
			intensityX += I*tmp[0]/10000/absr**3
			intensityY += I*tmp[1]/10000/absr**3
			intensityZ += I*tmp[2]/10000/absr**3
		pass
		return (intensityX,intensityY,intensityZ)
	pass
pass

class solenoid:
	# centerX is the center of X, R is the smallest radius of all solenoids
	def __init__(self,centerX,R,thickness,xNum,rNum,mesh,sign,I):
		self.I = I
		self.loops = []

		if (xNum%2) == 0:
			xNum = int(xNum/2)
			for j in range(xNum):
				for i in range(rNum):
					l = loop(centerX+(j+0.5)*thickness, R+i*thickness, sign, mesh)
					self.loops.append(l)
					l = loop(centerX-(j+0.5)*thickness, R+i*thickness, sign, mesh)
					self.loops.append(l)
		else:
			xNum = int((xNum+1)/2)
			for j in range(xNum):
				for i in range(rNum):
					if j == 0:
							l = loop(centerX, R+i*thickness, sign, mesh)
							self.loops.append(l)
					else:
						l = loop(centerX+j*thickness, R+i*thickness, sign, mesh)
						self.loops.append(l)
						l = loop(centerX-j*thickness, R+i*thickness, sign, mesh)
						self.loops.append(l)

	def addCoil(self,centerX,R,thickness,xNum,rNum,mesh,sign,I):
		if (xNum%2) == 0:
			xNum = int(xNum/2)
			for j in range(xNum):
				for i in range(rNum):
					l = loop(centerX+(j+0.5)*thickness, R+i*thickness, sign, mesh)
					self.loops.append(l)
					l = loop(centerX-(j+0.5)*thickness, R+i*thickness, sign, mesh)
					self.loops.append(l)
		else:
			xNum = int((xNum+1)/2)
			for j in range(xNum):
				for i in range(rNum):
					if j == 0:
							l = loop(centerX, R+i*thickness, sign, mesh)
							self.loops.append(l)
					else:
						l = loop(centerX+j*thickness, R+i*thickness, sign, mesh)
						self.loops.append(l)
						l = loop(centerX-j*thickness, R+i*thickness, sign, mesh)
						self.loops.append(l)


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
		surf = ax.plot_surface(x*100, z*100, inten*10, cmap=cm.coolwarm,
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
		surf = ax.plot_surface(y*100, z*100, inten*10, cmap=cm.coolwarm,
		                       linewidth=0, antialiased=False)
		fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()
		pass

	def gradientX(self,xMax):
		x = linspace(-xMax/100,xMax/100,100)
		intensity = self.SumField((x,0,0))
		for i in range(3):
			intensity[i] = gradient(intensity[i])
		plt.plot(x*100,intensity[0]*10,label="Bx")
		plt.plot(x*100,intensity[1]*10,label="By")
		plt.plot(x*100,intensity[2]*10,label="Bz")
		plt.xlabel('X(cm)')
		plt.ylabel('∂B/∂X(G/cm)')
		plt.title('Gradient of B vs position of X')
		plt.legend(loc=2)
		plt.show()
	def gradientZ(self,zMax):
		z = linspace(-zMax/100,zMax/100,100)
		intensity = self.SumField((0,0,z))
		for i in range(3):
			intensity[i] = gradient(intensity[i])
		plt.plot(z*100,intensity[0]*10,label="Bx")
		plt.plot(z*100,intensity[1]*10,label="By")
		plt.plot(z*100,intensity[2]*10,label="Bz")
		plt.xlabel('Z(cm)')
		plt.ylabel('∂B/∂Z(G/cm)')
		plt.title('Gradient of B vs position of Z')
		plt.legend(loc=2)
		plt.show()
	def BX(self,xMax):
		x = linspace(-xMax/100,xMax/100,100)
		intensity = self.SumField((x,0,0))
		plt.plot(x*100,intensity[0]*10,label="Bx")
		plt.plot(x*100,intensity[1]*10,label="By")
		plt.plot(x*100,intensity[2]*10,label="Bz")
		plt.xlabel('X(cm)')
		plt.ylabel('B(G)')
		plt.title('B vs position of X')
		plt.legend(loc=2)
		plt.show()
	def BZ(self,zMax):
		z = linspace(-zMax/100,zMax/100,100)
		intensity = self.SumField((0,0,z))
		plt.plot(z*100,intensity[0]*10,label="Bx")
		plt.plot(z*100,intensity[1]*10,label="By")
		plt.plot(z*100,intensity[2]*10,label="Bz")
		plt.xlabel('Z(cm)')
		plt.ylabel('B(G)')
		plt.title('B vs position of Z')
		plt.legend(loc=2)
		plt.show()


#centerX,R,thickness,xNum,rNum,mesh,sign,I
l = solenoid(-0.06, 0.01875, 0.0011,20 , 40, 50, 1,1.5)
l.addCoil(0.06, 0.01875, 0.0011, 20, 40, 50, -1, 1.5)
l.meshXZ(3, 5)
l.meshYZ(3, 3)
l.gradientX(10)
l.gradientZ(10)
l.BX(10)
l.BZ(10)
