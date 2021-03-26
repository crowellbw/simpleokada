#!/usr/bin/python
import math
import numpy

#This program computes the Green's functions from Okada's formulation for slip on fault patches at a set of station locations

eps = 6.1232e-14 #Small number to test cos(90) = 0
nu = 0.25 #Poisson's ratio

def greenF(e, n, depth, strike, dip, W, L):
	(l2,l1)=numpy.shape(e)
	G = numpy.zeros([3*l1,2*l2])
	for i in range (0, l2):
		for j in range (0, l1):
			strike1 = math.radians(strike[i])
			dip1 = math.radians(dip[i])
			d = depth[i,j] + math.sin(dip1)*W[i]/2
			ec = e[i,j] + math.cos(strike1)*math.cos(dip1)*W[i]/2
			nc = n[i,j] - math.sin(strike1)*math.cos(dip1)*W[i]/2
			x = math.cos(strike1)*nc + math.sin(strike1)*ec + L[i]/2
			y = math.sin(strike1)*nc - math.cos(strike1)*ec + math.cos(dip1)*W[i]
			
			p = y*math.cos(dip1) + d*math.sin(dip1);
			q = y*math.sin(dip1) - d*math.cos(dip1);

			g1 = -1/(2*math.pi) * (ux_ss(x,p,q,dip1) - ux_ss(x,p-W[i],q,dip1) - ux_ss(x-L[i],p,q,dip1) + ux_ss(x-L[i],p-W[i],q,dip1))
			g2 = -1/(2*math.pi) * (ux_ds(x,p,q,dip1) - ux_ds(x,p-W[i],q,dip1) - ux_ds(x-L[i],p,q,dip1) + ux_ds(x-L[i],p-W[i],q,dip1))

			g3 = -1/(2*math.pi) * (uy_ss(x,p,q,dip1) - uy_ss(x,p-W[i],q,dip1) - uy_ss(x-L[i],p,q,dip1) + uy_ss(x-L[i],p-W[i],q,dip1))
			g4 = -1/(2*math.pi) * (uy_ds(x,p,q,dip1) - uy_ds(x,p-W[i],q,dip1) - uy_ds(x-L[i],p,q,dip1) + uy_ds(x-L[i],p-W[i],q,dip1))

			g5 = -1/(2*math.pi) * (uz_ss(x,p,q,dip1) - uz_ss(x,p-W[i],q,dip1) - uz_ss(x-L[i],p,q,dip1) + uz_ss(x-L[i],p-W[i],q,dip1))
			g6 = -1/(2*math.pi) * (uz_ds(x,p,q,dip1) - uz_ds(x,p-W[i],q,dip1) - uz_ds(x-L[i],p,q,dip1) + uz_ds(x-L[i],p-W[i],q,dip1))


			g1n = math.sin(strike1)*g1 - math.cos(strike1)*g3
			g3n = math.cos(strike1)*g1 + math.sin(strike1)*g3

			g2n = math.sin(strike1)*g2 - math.cos(strike1)*g4
			g4n = math.cos(strike1)*g2 + math.sin(strike1)*g4

			G[3*j,2*i] = g1n
			G[3*j,2*i+1] = g2n

			G[3*j+1,2*i] = g3n
			G[3*j+1,2*i+1] = g4n

			G[3*j+2,2*i] = g5
			G[3*j+2,2*i+1] = g6
	return (G)		

def ux_ss(xi, eta, q, dip):
	R = math.sqrt(xi*xi+eta*eta+q*q)
	u = xi*q/(R*(R+eta)) + math.atan(xi*eta/(q*R)) + I1(xi, eta, q, dip, R)*math.sin(dip)
	
	return (u)

def uy_ss(xi, eta, q, dip):
	R = math.sqrt(xi*xi+eta*eta+q*q)
	yb = eta*math.cos(dip)+ q*math.sin(dip)
	u = yb*q/(R*(R+eta)) + q*math.cos(dip)/(R+eta) + I2(eta, q, dip, R)*math.sin(dip) 
	return (u)

def uz_ss(xi, eta, q, dip):
	R = math.sqrt(xi*xi+eta*eta+q*q)
	db = eta*math.sin(dip) - q*math.cos(dip)
	u = db*q/(R*(R+eta)) + q*math.sin(dip)/(R+eta) + I4(eta, q, dip, R)*math.sin(dip)

	return (u)

def ux_ds(xi, eta, q, dip):
	R = math.sqrt(xi*xi+eta*eta+q*q)
	u = q/R - I3(eta, q, dip, R)*math.sin(dip)*math.cos(dip)
	
	return (u)

def uy_ds(xi, eta, q, dip):
	R = math.sqrt(xi*xi+eta*eta+q*q)
	yb = eta*math.cos(dip)+ q*math.sin(dip)
	u = yb*q/(R*(R+xi)) + math.cos(dip)*math.atan(xi*eta/(q*R)) - I1(xi, eta, q, dip, R)*math.sin(dip)*math.cos(dip)
	
	return (u)

def uz_ds(xi, eta, q, dip):
	R = math.sqrt(xi*xi+eta*eta+q*q)
	db = eta*math.sin(dip) - q*math.cos(dip)
	u = db*q/(R*(R+xi)) + math.sin(dip)*math.atan(xi*eta/(q*R)) - I5(xi, eta, q, dip, R)*math.sin(dip)*math.cos(dip)
	
	return (u)

def I1(xi, eta, q, dip, R):
	db = eta*math.sin(dip) - q*math.cos(dip)
	if math.cos(dip) > eps:
		I = (1-2*nu)*(-xi/math.cos(dip)/(R+db)) - math.sin(dip)/math.cos(dip)*I5(xi, eta, q, dip, R)
		
	else:
		I = -(1-2*nu)/2*xi*q/math.pow(R+db,2)
	return (I)

def I2(eta, q, dip, R):
	I = (1-2*nu)*(-math.log(R+eta)) - I3(eta, q, dip, R)
	return (I)

def I3(eta, q, dip, R):
	yb = eta*math.cos(dip) + q*math.sin(dip)
	db = eta*math.sin(dip) - q*math.cos(dip)
	if math.cos(dip) > eps:
		I = (1-2*nu)*(yb/math.cos(dip)/(R+db)-math.log(R+eta)) + math.sin(dip)/math.cos(dip)*I4(eta, q, dip, R)
	else:
		I = (1-2*nu)/2*( eta/(R+db) + yb*q/math.pow(R+db,2) - math.log(R+eta) )
	return (I)

def I4(eta, q, dip, R):
	db = eta*math.sin(dip) - q*math.cos(dip)
	if math.cos(dip) > eps:
		I = (1-2*nu)/math.cos(dip) * (math.log(R+db) - math.sin(dip)*math.log(R+eta))
	else:
		I = -(1-2*nu)*q/(R+db)
	return (I)

def I5(xi, eta, q, dip, R):
	db = eta*math.sin(dip) - q*math.cos(dip)
	X = math.sqrt(xi*xi+q*q)
	if math.cos(dip) > eps:
		I = (1-2*nu)*2/math.cos(dip)* math.atan( ( eta*(X+q*math.cos(dip)) + X*(R+X)*math.sin(dip) )/(xi*(R+X)*math.cos(dip)))
	else:
		I = -(1-2*nu)*xi*math.sin(dip)/(R+db)
	return (I)





