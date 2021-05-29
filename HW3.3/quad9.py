from numpy import array, sqrt, zeros, ix_
from scipy.linalg import det, inv


def quad9(xy, properties):

	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	x3 = xy[3,0]
	x4 = xy[4,0]
	x5 = xy[5,0]
	x6 = xy[6,0]
	x7 = xy[7,0]
	x8 = xy[8,0]
    
	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]
	y3 = xy[3,1]
	y4 = xy[4,1]
	y5 = xy[5,1]
	y6 = xy[6,1]
	y7 = xy[7,1]
	y8 = xy[8,1]

	ke = zeros((18,18))
	fe = zeros((18,1))

	#Primer punto de Gauss de la regla 2x2
	# xi = 1.0 / sqrt(3)
	# eta = -1.0 / sqrt(3)
	# wi = 1.0
	# wj = 1.0

	gauss_rule = [
		( -sqrt(3/5), -sqrt(3/5), 5/9, 5/9),               #0
		(  sqrt(3/5), -sqrt(3/5), 5/9, 5/9),               #1
		(  sqrt(3/5),  sqrt(3/5), 5/9, 5/9),               #2
		(- sqrt(3/5),  sqrt(3/5), 5/9, 5/9),               #3
		(  0.0      , -sqrt(3/5), 8/9, 5/9),               #4*
		(  sqrt(3/5), 0.0       , 5/9, 8/9),               #5
		( 0.0       ,  sqrt(3/5), 8/9, 5/9),               #6
		( -sqrt(3/5), 0.0       , 5/9, 8/9),               #7
		(  0.0      ,     0.0   , 8/9, 8/9)                #8
	]

	for ξ, η, wi, wj in gauss_rule:

		#print(f"ξ = {ξ} η = {η}")

		x = x0*η*ξ*(η - 1)*(ξ - 1)/4 + x1*η*ξ*(η - 1)*(ξ + 1)/4 + x2*η*ξ*(η + 1)*(ξ + 1)/4 + x3*η*ξ*(η + 1)*(ξ - 1)/4 + x4*η*(1 - ξ**2)*(η - 1)/2 + x5*ξ*(1 - η**2)*(ξ + 1)/2 + x6*η*(1 - ξ**2)*(η + 1)/2 + x7*ξ*(1 - η**2)*(ξ - 1)/2 + x8*(1 - η**2)*(1 - ξ**2)
		y = y0*η*ξ*(η - 1)*(ξ - 1)/4 + y1*η*ξ*(η - 1)*(ξ + 1)/4 + y2*η*ξ*(η + 1)*(ξ + 1)/4 + y3*η*ξ*(η + 1)*(ξ - 1)/4 + y4*η*(1 - ξ**2)*(η - 1)/2 + y5*ξ*(1 - η**2)*(ξ + 1)/2 + y6*η*(1 - ξ**2)*(η + 1)/2 + y7*ξ*(1 - η**2)*(ξ - 1)/2 + y8*(1 - η**2)*(1 - ξ**2)
		dx_dxi = x0*η*ξ*(η - 1)/4 + x0*η*(η - 1)*(ξ - 1)/4 + x1*η*ξ*(η - 1)/4 + x1*η*(η - 1)*(ξ + 1)/4 + x2*η*ξ*(η + 1)/4 + x2*η*(η + 1)*(ξ + 1)/4 + x3*η*ξ*(η + 1)/4 + x3*η*(η + 1)*(ξ - 1)/4 - x4*η*ξ*(η - 1) + x5*ξ*(1 - η**2)/2 + x5*(1 - η**2)*(ξ + 1)/2 - x6*η*ξ*(η + 1) + x7*ξ*(1 - η**2)/2 + x7*(1 - η**2)*(ξ - 1)/2 - 2*x8*ξ*(1 - η**2)
		dx_deta = x0*η*ξ*(ξ - 1)/4 + x0*ξ*(η - 1)*(ξ - 1)/4 + x1*η*ξ*(ξ + 1)/4 + x1*ξ*(η - 1)*(ξ + 1)/4 + x2*η*ξ*(ξ + 1)/4 + x2*ξ*(η + 1)*(ξ + 1)/4 + x3*η*ξ*(ξ - 1)/4 + x3*ξ*(η + 1)*(ξ - 1)/4 + x4*η*(1 - ξ**2)/2 + x4*(1 - ξ**2)*(η - 1)/2 - x5*η*ξ*(ξ + 1) + x6*η*(1 - ξ**2)/2 + x6*(1 - ξ**2)*(η + 1)/2 - x7*η*ξ*(ξ - 1) - 2*x8*η*(1 - ξ**2)
		dy_dxi = y0*η*ξ*(η - 1)/4 + y0*η*(η - 1)*(ξ - 1)/4 + y1*η*ξ*(η - 1)/4 + y1*η*(η - 1)*(ξ + 1)/4 + y2*η*ξ*(η + 1)/4 + y2*η*(η + 1)*(ξ + 1)/4 + y3*η*ξ*(η + 1)/4 + y3*η*(η + 1)*(ξ - 1)/4 - y4*η*ξ*(η - 1) + y5*ξ*(1 - η**2)/2 + y5*(1 - η**2)*(ξ + 1)/2 - y6*η*ξ*(η + 1) + y7*ξ*(1 - η**2)/2 + y7*(1 - η**2)*(ξ - 1)/2 - 2*y8*ξ*(1 - η**2)
		dy_deta = y0*η*ξ*(ξ - 1)/4 + y0*ξ*(η - 1)*(ξ - 1)/4 + y1*η*ξ*(ξ + 1)/4 + y1*ξ*(η - 1)*(ξ + 1)/4 + y2*η*ξ*(ξ + 1)/4 + y2*ξ*(η + 1)*(ξ + 1)/4 + y3*η*ξ*(ξ - 1)/4 + y3*ξ*(η + 1)*(ξ - 1)/4 + y4*η*(1 - ξ**2)/2 + y4*(1 - ξ**2)*(η - 1)/2 - y5*η*ξ*(ξ + 1) + y6*η*(1 - ξ**2)/2 + y6*(1 - ξ**2)*(η + 1)/2 - y7*η*ξ*(ξ - 1) - 2*y8*η*(1 - ξ**2)

		dN0_dxi  =  η*ξ*(η - 1)/4 + η*(η - 1)*(ξ - 1)/4
		dN0_deta =  η*ξ*(ξ - 1)/4 + ξ*(η - 1)*(ξ - 1)/4
		dN1_dxi  =  η*ξ*(η - 1)/4 + η*(η - 1)*(ξ + 1)/4
		dN1_deta =  η*ξ*(ξ + 1)/4 + ξ*(η - 1)*(ξ + 1)/4
		dN2_dxi  =  η*ξ*(η + 1)/4 + η*(η + 1)*(ξ + 1)/4
		dN2_deta =  η*ξ*(ξ + 1)/4 + ξ*(η + 1)*(ξ + 1)/4
		dN3_dxi  =  η*ξ*(η + 1)/4 + η*(η + 1)*(ξ - 1)/4
		dN3_deta =  η*ξ*(ξ - 1)/4 + ξ*(η + 1)*(ξ - 1)/4
		dN4_dxi  = -η*ξ*(η - 1)
		dN4_deta =  η*(1 - ξ**2)/2 + (1/2 - ξ**2/2)*(η - 1)
		dN5_dxi  =  ξ*(1 - η**2)/2 + (1 - η**2)*(ξ/2 + 1/2)
		dN5_deta = -η*ξ*(ξ + 1)
		dN6_dxi  = -η*ξ*(η + 1)
		dN6_deta =  η*(1 - ξ**2)/2 + (1 - ξ**2)*(η/2 + 1/2)
		dN7_dxi  =  ξ*(1 - η**2)/2 + (1/2 - η**2/2)*(ξ - 1)
		dN7_deta = -η*ξ*(ξ - 1)
		dN8_dxi  = -2*ξ*(1/2 - η**2/2)
		dN8_deta = -η*(1 - ξ**2)*2

		#print(f"x = {x} y = {y}")


		J = array([
		[dx_dxi, dx_deta],
		[dy_dxi, dy_deta]
		]).T

		detJ = det(J)
		#print(f"det(J)={detJ}")
		if detJ <= 0.:
			print("FATAL! detJ <= 0...")
			exit(-1)

		Jinv = inv(J)

		# print(f"J = {J}")
		# print(f"detJ = {detJ}")

		dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
		dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
		dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
		dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])
		dN4_dxy = Jinv@array([ dN4_dxi, dN4_deta ])
		dN5_dxy = Jinv@array([ dN5_dxi, dN5_deta ])
		dN6_dxy = Jinv@array([ dN6_dxi, dN6_deta ])
		dN7_dxy = Jinv@array([ dN7_dxi, dN7_deta ])
		dN8_dxy = Jinv@array([ dN8_dxi, dN8_deta ])


		# ε = B ue
		B = zeros((3, 18))

		B[0,0] = dN0_dxy[0]
		B[1,1] = dN0_dxy[1]
		B[2,0] = dN0_dxy[1]
		B[2,1] = dN0_dxy[0]
        
		B[0,2] = dN1_dxy[0]
		B[1,3] = dN1_dxy[1]
		B[2,2] = dN1_dxy[1]
		B[2,3] = dN1_dxy[0]
        
		B[0,4] = dN2_dxy[0]
		B[1,5] = dN2_dxy[1]
		B[2,4] = dN2_dxy[1]
		B[2,5] = dN2_dxy[0]
        
		B[0,6] = dN3_dxy[0]
		B[1,7] = dN3_dxy[1]
		B[2,6] = dN3_dxy[1]
		B[2,7] = dN3_dxy[0]
        
		B[0,8] = dN4_dxy[0]
		B[1,9] = dN4_dxy[1]
		B[2,8] = dN4_dxy[1]
		B[2,9] = dN4_dxy[0]
        
		B[0,10] = dN5_dxy[0]
		B[1,11] = dN5_dxy[1]
		B[2,10] = dN5_dxy[1]
		B[2,11] = dN5_dxy[0]
        
		B[0,12] = dN6_dxy[0]
		B[1,13] = dN6_dxy[1]
		B[2,12] = dN6_dxy[1]
		B[2,13] = dN6_dxy[0]
        
		B[0,14] = dN7_dxy[0]
		B[1,15] = dN7_dxy[1]
		B[2,14] = dN7_dxy[1]
		B[2,15] = dN7_dxy[0]
        
		B[0,16] = dN8_dxy[0]
		B[1,17] = dN8_dxy[1]
		B[2,16] = dN8_dxy[1]
		B[2,17] = dN8_dxy[0]
        
		#print(f"B = {B}")


		ke += t * wi * wj * B.T @ Eσ @ B * detJ

		#print(f"ke = {ke}")
	return ke, fe
    
def quad9_line_load(xy, properties_load):
    
    e = properties_load["t"]
    tx = properties_load["tx"]
    ty = properties_load["ty"]
    ke = 0
    
    x0 = xy[0,0]
    x1 = xy[1,0]
    x2 = xy[2,0]
    y0 = xy[0,1]
    y1 = xy[1,1]
    y2 = xy[2,1]
    L1 = sqrt((x2 - x0)**2 + (y2 - y0)**2)
    L2 = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    t1 = array([tx,ty,0,0,tx,ty])
    t2 = array([0,0,tx,ty,tx,ty])
    fe = e*L1/2*t1 + e*L2/2*t2
    print(fe)
    return ke, fe
    
def quad9_post(xy, u_e, properties):

	E = properties["E"]
	ν = properties["nu"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	#Podemos pasarle otros valores de xi y eta
	if "xi" in properties:
		ξ = properties["xi"]
	else:
		ξ = 0.0                        ##?

	if "eta" in properties:
		η = properties["eta"]
	else:
		η = 0.0                        ##?

	Eσ = E / (1-ν**2) * array(
		[
		[1 , ν , 0       ]       ,
		[ν , 1 , 0       ]       ,
		[0 , 0 , (1-ν)/2 ]
		])

	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	x3 = xy[3,0]
	x4 = xy[4,0]
	x5 = xy[5,0]
	x6 = xy[6,0]
	x7 = xy[7,0]
	x8 = xy[8,0]
    
	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]
	y3 = xy[3,1]
	y4 = xy[4,1]
	y5 = xy[5,1]
	y6 = xy[6,1]
	y7 = xy[7,1]
	y8 = xy[8,1]

	x = x0*η*ξ*(η - 1)*(ξ - 1)/4 + x1*η*ξ*(η - 1)*(ξ + 1)/4 + x2*η*ξ*(η + 1)*(ξ + 1)/4 + x3*η*ξ*(η + 1)*(ξ - 1)/4 + x4*η*(1 - ξ**2)*(η - 1)/2 + x5*ξ*(1 - η**2)*(ξ + 1)/2 + x6*η*(1 - ξ**2)*(η + 1)/2 + x7*ξ*(1 - η**2)*(ξ - 1)/2 + x8*(1 - η**2)*(1 - ξ**2)
	y = y0*η*ξ*(η - 1)*(ξ - 1)/4 + y1*η*ξ*(η - 1)*(ξ + 1)/4 + y2*η*ξ*(η + 1)*(ξ + 1)/4 + y3*η*ξ*(η + 1)*(ξ - 1)/4 + y4*η*(1 - ξ**2)*(η - 1)/2 + y5*ξ*(1 - η**2)*(ξ + 1)/2 + y6*η*(1 - ξ**2)*(η + 1)/2 + y7*ξ*(1 - η**2)*(ξ - 1)/2 + y8*(1 - η**2)*(1 - ξ**2)
	dx_dxi = x0*η*ξ*(η - 1)/4 + x0*η*(η - 1)*(ξ - 1)/4 + x1*η*ξ*(η - 1)/4 + x1*η*(η - 1)*(ξ + 1)/4 + x2*η*ξ*(η + 1)/4 + x2*η*(η + 1)*(ξ + 1)/4 + x3*η*ξ*(η + 1)/4 + x3*η*(η + 1)*(ξ - 1)/4 - x4*η*ξ*(η - 1) + x5*ξ*(1 - η**2)/2 + x5*(1 - η**2)*(ξ + 1)/2 - x6*η*ξ*(η + 1) + x7*ξ*(1 - η**2)/2 + x7*(1 - η**2)*(ξ - 1)/2 - 2*x8*ξ*(1 - η**2)
	dx_deta = x0*η*ξ*(ξ - 1)/4 + x0*ξ*(η - 1)*(ξ - 1)/4 + x1*η*ξ*(ξ + 1)/4 + x1*ξ*(η - 1)*(ξ + 1)/4 + x2*η*ξ*(ξ + 1)/4 + x2*ξ*(η + 1)*(ξ + 1)/4 + x3*η*ξ*(ξ - 1)/4 + x3*ξ*(η + 1)*(ξ - 1)/4 + x4*η*(1 - ξ**2)/2 + x4*(1 - ξ**2)*(η - 1)/2 - x5*η*ξ*(ξ + 1) + x6*η*(1 - ξ**2)/2 + x6*(1 - ξ**2)*(η + 1)/2 - x7*η*ξ*(ξ - 1) - 2*x8*η*(1 - ξ**2)
	dy_dxi = y0*η*ξ*(η - 1)/4 + y0*η*(η - 1)*(ξ - 1)/4 + y1*η*ξ*(η - 1)/4 + y1*η*(η - 1)*(ξ + 1)/4 + y2*η*ξ*(η + 1)/4 + y2*η*(η + 1)*(ξ + 1)/4 + y3*η*ξ*(η + 1)/4 + y3*η*(η + 1)*(ξ - 1)/4 - y4*η*ξ*(η - 1) + y5*ξ*(1 - η**2)/2 + y5*(1 - η**2)*(ξ + 1)/2 - y6*η*ξ*(η + 1) + y7*ξ*(1 - η**2)/2 + y7*(1 - η**2)*(ξ - 1)/2 - 2*y8*ξ*(1 - η**2)
	dy_deta = y0*η*ξ*(ξ - 1)/4 + y0*ξ*(η - 1)*(ξ - 1)/4 + y1*η*ξ*(ξ + 1)/4 + y1*ξ*(η - 1)*(ξ + 1)/4 + y2*η*ξ*(ξ + 1)/4 + y2*ξ*(η + 1)*(ξ + 1)/4 + y3*η*ξ*(ξ - 1)/4 + y3*ξ*(η + 1)*(ξ - 1)/4 + y4*η*(1 - ξ**2)/2 + y4*(1 - ξ**2)*(η - 1)/2 - y5*η*ξ*(ξ + 1) + y6*η*(1 - ξ**2)/2 + y6*(1 - ξ**2)*(η + 1)/2 - y7*η*ξ*(ξ - 1) - 2*y8*η*(1 - ξ**2)

	dN0_dxi  =  η*ξ*(η - 1)/4 + η*(η - 1)*(ξ - 1)/4
	dN0_deta =  η*ξ*(ξ - 1)/4 + ξ*(η - 1)*(ξ - 1)/4
	dN1_dxi  =  η*ξ*(η - 1)/4 + η*(η - 1)*(ξ + 1)/4
	dN1_deta =  η*ξ*(ξ + 1)/4 + ξ*(η - 1)*(ξ + 1)/4
	dN2_dxi  =  η*ξ*(η + 1)/4 + η*(η + 1)*(ξ + 1)/4
	dN2_deta =  η*ξ*(ξ + 1)/4 + ξ*(η + 1)*(ξ + 1)/4
	dN3_dxi  =  η*ξ*(η + 1)/4 + η*(η + 1)*(ξ - 1)/4
	dN3_deta =  η*ξ*(ξ - 1)/4 + ξ*(η + 1)*(ξ - 1)/4
	dN4_dxi  = -η*ξ*(η - 1)
	dN4_deta =  η*(1 - ξ**2)/2 + (1/2 - ξ**2/2)*(η - 1)
	dN5_dxi  =  ξ*(1 - η**2)/2 + (1 - η**2)*(ξ/2 + 1/2)
	dN5_deta = -η*ξ*(ξ + 1)
	dN6_dxi  = -η*ξ*(η + 1)
	dN6_deta =  η*(1 - ξ**2)/2 + (1 - ξ**2)*(η/2 + 1/2)
	dN7_dxi  =  ξ*(1 - η**2)/2 + (1/2 - η**2/2)*(ξ - 1)
	dN7_deta = -η*ξ*(ξ - 1)
	dN8_dxi  = -2*ξ*(1/2 - η**2/2)
	dN8_deta = -η*(1 - ξ**2)*2
	# print(f"x = {x} y = {y}")

	# print(f"x = {x} y = {y}")

	J = array([
	[dx_dxi, dx_deta],
	[dy_dxi, dy_deta]
	]).T

	detJ = det(J)

	if detJ <= 0.:
		print(f"FATAL! detJ <= 0...")
		exit(-1)

	Jinv = inv(J)

	# print(f"J = {J}")
	# print(f"detJ = {detJ}")

	dN0_dxy = Jinv@array([ dN0_dxi, dN0_deta ])
	dN1_dxy = Jinv@array([ dN1_dxi, dN1_deta ])
	dN2_dxy = Jinv@array([ dN2_dxi, dN2_deta ])
	dN3_dxy = Jinv@array([ dN3_dxi, dN3_deta ])
	dN4_dxy = Jinv@array([ dN4_dxi, dN4_deta ])
	dN5_dxy = Jinv@array([ dN5_dxi, dN5_deta ])
	dN6_dxy = Jinv@array([ dN6_dxi, dN6_deta ])
	dN7_dxy = Jinv@array([ dN7_dxi, dN7_deta ])
	dN8_dxy = Jinv@array([ dN8_dxi, dN8_deta ])

	# ε = B ue
	B = zeros((3, 18))
	B[0,0] = dN0_dxy[0]
	B[1,1] = dN0_dxy[1]
	B[2,0] = dN0_dxy[1]
	B[2,1] = dN0_dxy[0]
	B[0,2] = dN1_dxy[0]
	B[1,3] = dN1_dxy[1]
	B[2,2] = dN1_dxy[1]
	B[2,3] = dN1_dxy[0]
	B[0,4] = dN2_dxy[0]
	B[1,5] = dN2_dxy[1]
	B[2,4] = dN2_dxy[1]
	B[2,5] = dN2_dxy[0]
	B[0,6] = dN3_dxy[0]
	B[1,7] = dN3_dxy[1]
	B[2,6] = dN3_dxy[1]
	B[2,7] = dN3_dxy[0]
       
	B[0,8] = dN4_dxy[0]
	B[1,9] = dN4_dxy[1]
	B[2,8] = dN4_dxy[1]
	B[2,9] = dN4_dxy[0]
        
	B[0,10] = dN5_dxy[0]
	B[1,11] = dN5_dxy[1]
	B[2,10] = dN5_dxy[1]
	B[2,11] = dN5_dxy[0]
        
	B[0,12] = dN6_dxy[0]
	B[1,13] = dN6_dxy[1]
	B[2,12] = dN6_dxy[1]
	B[2,13] = dN6_dxy[0]
       
	B[0,14] = dN7_dxy[0]
	B[1,15] = dN7_dxy[1]
	B[2,14] = dN7_dxy[1]
	B[2,15] = dN7_dxy[0]
       
	B[0,16] = dN8_dxy[0]
	B[1,17] = dN8_dxy[1]
	B[2,16] = dN8_dxy[1]
	B[2,17] = dN8_dxy[0]


	ε = B @ u_e
	σ = Eσ @ ε

	return ε, σ

# xy = array([
# [-1,-1],
# [1,-1],
# [1,1],
# [-1,1],
# 	])

# xy = array([
# [0,0],
# [1,0],
# [1,1],
# [0,1],
# 	])

# properties = {}
# properties["E"] = 1.
# properties["nu"] = 0.25
# properties["bx"] = 0
# properties["by"] = 1.
# properties["t"] = 1.

# ke, fe = quad4(xy, properties)

# print(f"ke = {ke}")

# fixed_dofs = [0, 1, 2, 3]
# free_dofs = [4, 5, 6, 7]

# ke_ff = ke[ix_(free_dofs, free_dofs)]
# fe_ff = array([0, -1, 0, -1])

# print(f"ke_ff = {ke_ff}")

# from scipy.linalg import solve

# u = zeros((8,1))
# uf = solve(ke_ff, fe_ff)

# print(f"uf = {uf}")