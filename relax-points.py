# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 05:56:51 2024

It creates an even sampling of points in the surface of a sphere. 
Starts with a Fibonacci distribution of points and proceed 
with a Spherical Voronoi relaxation (see docs)
It opens a folder named with the number of samples and save the
converged points and the last Voronoi diagram.

@author: Nicolas Passarelli
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import SphericalVoronoi, geometric_slerp
import os


		############################################################## input
samples=20 #numer of points 
nitersmax=10000 # max iterations
step=0.2  #for the voronoi relaxation
tol=0.0005 # convergence for the RSTD
steptol=10 # how many iterations are considered in RSTD

radius0 = 1 #radius of big sphere
center = np.array([0, 0, 0]) # center of sphere

plotFibonacci=True #initial sampling plot
dovideo=True #save frames for making video

		############################################################## main
#make directory for output
dirn='samples='+ str(samples)
dirv='vid'
if not os.path.exists(dirn):
    os.mkdir(dirn)
	
if not os.path.exists(dirv):
    os.mkdir(dirv)


dirn=dirn +'/' 
		
delta=1./float(samples-1)
phi=2/(np.sqrt(5)+1)

def fibonacci_sphere(samples):
	points=[]		
	for i in range(samples):
		y=1-delta*2*i
		radius=np.sqrt(radius0**2-y**2)
		theta=phi*i
		x=np.cos(theta)*radius
		z=np.sin(theta)*radius
		points.append([x,y,z])
	return points

def fibonacci_sphere2(samples):
	points=[]
	for i in range(samples):
		z=1-delta*2*i
		radius=np.sqrt(radius0**2-z**2)
		theta=phi*i
		x=np.cos(theta)*radius
		y=np.sin(theta)*radius
		points.append([x,y,z])
	return points


def length(a): return np.sqrt(np.dot(a,a))

def angleBetweenUnitVectors(a,b):
	# https://www.plunk.org/~hatch/rightway.html
	if np.dot(a,b) < 0.:
		return np.pi - 2*np.arcsin(length(-b-a)/2.)
	else:
		return 2*np.arcsin(length(b-a)/2.)
	
		##############################################################init
centers=fibonacci_sphere2(samples)
centers=np.asarray(centers)	
ncenters=centers.shape[0]
rad=0.1 #rad for plot

xs=centers[:,0]
ys=centers[:,1]
zs=centers[:,2]

rs=np.sqrt(centers[:,0]**2+centers[:,1]**2+centers[:,2]**2)
thetas=np.arctan2(ys,xs)
phis=np.arccos(zs/rs)

if plotFibonacci:
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	ax.scatter(xs, ys, zs)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_xlim([-1,1])
	ax.set_ylim([-1,1])
	ax.set_zlim([-1,1])
	ax.set_aspect('equal')
	plt.show()

	fig,ax = plt.subplots()
	ax.scatter(thetas, phis)
	ax.set_xlabel('$\\theta$')
	ax.set_ylabel('$\\phi$')
	ax.set_xlim([-np.pi,np.pi])
	ax.set_ylim([0,np.pi])
	ax.set_aspect('equal')
	fig.set_size_inches(4, 3)
	plt.tight_layout()
	plt.savefig('./'+dirn+'Fibonacci2.png')

	plt.show()


	######################################################### loop
	
points=centers #points moves in optimisation

frames=[]
listRSD=[]
jj=0
while jj<nitersmax:
	jj+=1
	
	######################################################### voronoi
	sv = SphericalVoronoi(points, radius0, center, threshold=1e-06)
	areas=sv.calculate_areas()
	areasRSTD=np.std(areas)/np.mean(areas)
	
	
	centroid=np.zeros((len(sv.regions),3))
	for ii,region in enumerate(sv.regions):
		nv = len(region)	
		C=[0,0,0]
		for i in range(nv):
			v1=	sv.vertices[region[i]]
			v2=	sv.vertices[region[(i + 1) % nv]]
			v1u= v1/length(v1) 
			v2u= v2/length(v2)  		
			Na= np.cross(v1u,v2u)	# normal to the big circle, plane that pases throug origin		
			Nau= Na/ length(Na)
			dab=angleBetweenUnitVectors(v2u,v1u)
			C+=Nau*dab/2
		centroid[ii,0:3]= C / length(C) *radius0 
		angletofirst = angleBetweenUnitVectors(points[ii],centroid[ii,:])
		if (angletofirst > np.pi/2): centroid[ii,:] = - centroid[ii,:]
			
			
			#################################################################plot
	
	# sort vertices (optional, helpful for plotting)
	sv.sort_vertices_of_regions()
	t_vals = np.linspace(0, 1, 2000)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	# plot the unit sphere for reference (optional)
	u = np.linspace(0, 2 * np.pi, 100)
	v = np.linspace(0, np.pi, 100)
	x = np.outer(np.cos(u), np.sin(v))
	y = np.outer(np.sin(u), np.sin(v))
	z = np.outer(np.ones(np.size(u)), np.cos(v))
	ax.plot_surface(x, y, z, color='y', alpha=.3)
	# plot generator points

	ax.azim = 0
	ax.elev = 45
	_ = ax.set_xticks([])
	_ = ax.set_yticks([])
	_ = ax.set_zticks([])
	ax.set_aspect('equal')
	fig.set_size_inches(4, 4)
		
	
	ff=ax.scatter(centroid[:, 0], centroid[:, 1], centroid[:, 2], c='magenta',marker='x')
	ff=ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2],s=1, c='g')
	ff=ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')

	for region in sv.regions:
		n = len(region)
		for i in range(n):
			start = sv.vertices[region][i]
			end = sv.vertices[region][(i + 1) % n]
			result = geometric_slerp(start, end, t_vals)
			ff=ax.plot(result[..., 0],
					result[..., 1],
					result[..., 2],
					c='g',lw=0.1)
	
	stri='./vid/'+str(jj)+'.png'
	if dovideo: plt.savefig(stri)
	if dovideo: frames.append(stri)
	plt.show()
		
		
			################################################################# calculate next step
	points2=points*1
	for i in range(ncenters):
		points2[i,:]=geometric_slerp(points[i,:], centroid[i,:], step) 
	points=points2*1
	listRSD.append(areasRSTD)
	if jj>steptol:			#finish while loop when converges		
	
		last=listRSD[-1:-steptol:-1]	
		rstd_rstd=np.std(last)/np.mean(last)
		print(jj, 'area convergence',areasRSTD,'/',tol)
		print(jj, 'sigma convergence',rstd_rstd,'/',tol)
		
		if rstd_rstd<tol:
			jj=nitersmax
			
		if areasRSTD<tol:
			jj=nitersmax 		
			
			################################################################# save
			
np.savetxt('./'+dirn+'centers.dat',points2)
np.savetxt('./'+dirn+'vertices.dat',sv.vertices)


with open('./'+dirn+'faces.dat', 'w') as f:
	for region in sv.regions:
		st="\t".join( map( str, region ) ) + '\n'
		f.write(st)
	

if dovideo:
	import subprocess

	
	with open('.makevideo.sh', 'w') as f:
		f.write("magick ./vid/\%d.png[1-"+ str(len(frames)) +"] mymovie.gif")
	subprocess.call(['.\.makevideo.sh'], shell=True)
	
	

	 
	
	
	


			