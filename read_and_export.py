# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 05:56:51 2024

Read the results saved in folder.
It takes the voronoi diagram and create obj1 and save a .obj
It calculates the convex hull and create obj2 and save a .obj
see docs

@author: Nicolas Passarelli
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import  geometric_slerp, ConvexHull
  
#################################################################input

samples=20  #gives the dirname and number of faces of 
plotAngDist=True #angular plot of results

#################################################################plot


dirn='samples='+ str(samples)+'/'
points=np.genfromtxt('./'+dirn+'centers.dat',)
vertices=np.genfromtxt('./'+dirn+'vertices.dat')

#read faces from our own format
faces=[]
with open('./'+dirn+'faces.dat', 'r') as f:
	for line in f.readlines():
		sp=line.split(sep='\t')
		lis=[]
		for s in sp:
			lis.append(int(s))
		faces.append(lis)	
	
		

t_vals = np.linspace(0, 1, 2000)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='y', alpha=.6)

for region in faces:
	n = len(region)
	for i in range(n):
		start = vertices[region][i]
		end = vertices[region][(i + 1) % n]
		result = geometric_slerp(start, end, t_vals)
		ax.plot(result[..., 0],
				result[..., 1],
				result[..., 2],
				c='k',lw=0.5)

ax.azim = 60
ax.elev = 60
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_aspect('equal')
fig.set_size_inches(4, 4)
plt.show()


if plotAngDist:
	R=(points[:,0]**2+points[:,1]**2+points[:,2]**2)**.5 
	thetas=np.arctan2(points[:,1],points[:,0])
	phis=np.arccos(points[:,2]/R)

	R2=(vertices[:,0]**2+vertices[:,1]**2+vertices[:,2]**2)**.5 
	thetas2=np.arctan2(vertices[:,1],vertices[:,0])
	phis2=np.arccos(vertices[:,2]/R2)
	verticesA=np.vstack([thetas2,phis2])
	
	fig,ax = plt.subplots()
	ax.scatter(thetas, phis,label='vtx Obj 2')
	ax.scatter(thetas2, phis2,label='vtx Obj 1')
	
	for region in faces:
		n = len(region)
		for i in range(n):
			start = vertices[region][i,:]
			end = vertices[region][(i + 1) % n,:]		
			
			t1=np.arctan2(start[1],start[0])
			t2=np.arctan2(end[1],end[0])
			r1=(start[0]**2+start[1]**2+start[2]**2)**.5 
			r2=(end[0]**2+end[1]**2+end[2]**2)**.5 
			p1=np.arccos(start[2]/r1)
			p2=np.arccos(end[2]/r2)	
			
			#periodic boundaries handling
			if t2-t1 > np.pi :
				plt.plot((t1+2*np.pi,t2), (p1,p2), 'k', lw=0.5)	
				t2-=2*np.pi				
				
			if t2-t1 <-np.pi :
				plt.plot((t1-2*np.pi,t2), (p1,p2), 'k', lw=0.5)	
				t2+=2*np.pi 	

			if p2-p1 > np.pi/2 :
				plt.plot((t1,t2), (p1+2*np.pi,p2), 'k', lw=0.5)	
				p2-=2*np.pi 
				
			if p2-p1 <-np.pi/2 :
				plt.plot((t1,t2), (p1-2*np.pi,p2), 'k', lw=0.5)	
				p2+=2*np.pi 				

			plt.plot((t1,t2), (p1,p2), 'k', lw=0.5)	
			
	ax.set_xlabel('$\\theta$')
	ax.set_ylabel('$\\phi$')
	ax.set_xlim([-np.pi,np.pi])
	ax.set_ylim([0,np.pi])
	#plt.legend(loc='top left')
	ax.set_aspect('equal')
	fig.set_size_inches(4, 3)
	plt.tight_layout()
	plt.savefig('./'+dirn+'Angles.png')
	plt.show()	
	
	
#save obj file

def saveobj(filen, vertices, faces):
		
	with open(filen, 'w') as f:
		s1='# List of geometric vertices, with (x, y, z, [w]) coordinates, w is optional and defaults to 1.0.'
		f.write(s1)
		f.write('\n')
		
		for i in range(len(vertices)):
			f.write('v ')
			for j in range(3):	
				f.write(str(vertices[i,j])+' ')
			f.write('\n')
			
		s4='# Polygonal face element'
		f.write(s4)
		f.write('\n')		
		for region in faces:
			f.write('f ')
			for v in region:	
				f.write(str(v+1)+' ' ) 
			f.write(' \n')		
		s5='# Line element'
		f.write(s5)
		f.write('\n')
		for region in faces:
			f.write('l ')
			for v in region:	
				f.write(str(v+1)+' ' )
			f.write(' \n')	
		
#save	
filen='./'+dirn+'objfile1_'+str(samples) +'faces' +'.obj'
saveobj(filen, vertices, faces)

	
################################################################# convex Hull	
		
ch=ConvexHull(points)

vertices=ch.points
faces=ch.simplices

filen='./'+dirn+'objfile2_'+str(len(faces)) +'faces' +'.obj'
saveobj(filen, vertices, faces)


#################################################################plot

t_vals = np.linspace(0, 1, 2000)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='y', alpha=.6)


for region in faces:
	n = len(region)
	for i in range(n):
		start = vertices[region][i]
		end = vertices[region][(i + 1) % n]
		result = geometric_slerp(start, end, t_vals)
		ax.plot(result[..., 0],
				result[..., 1],
				result[..., 2],
				c='k',lw=0.5)

ax.azim = 60
ax.elev = 60
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_aspect('equal')
fig.set_size_inches(4, 4)
plt.show()

