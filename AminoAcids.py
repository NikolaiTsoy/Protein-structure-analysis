from __future__ import print_function
from pymol import cmd

import numpy as np
from math import acos, pi, sqrt





def get_vec_point(model, chains, res1, res2):
	chainLet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'O']
	
	cmd.fetch(model, async=0)

	#cmd.copy('{}_save'.format(model), model)


	#Initializ arguments
	Mean_point = [0,0,0]

	p1Sum = 0
	p2Sum = 0
	p3Sum = 0


	angM = 0
	sumErr = 0

	e1MeanSum = 0
	e2MeanSum = 0
	e3MeanSum = 0

	Mean_axis = [0,0,0]
	Mean_point = [0,0,0]

	#Calcul of mean point, vector with loop of aminoacid residues

	for q in range(int(res2)-int(res1)):
		res = int(res1) + q

		a = cmd.get_coords('/{}//{}/{}/CA'.format(model, chainLet[0], str(res)))[0]

		b = cmd.get_coords('/{}//{}/{}/CA'.format(model, chainLet[1], str(res)))[0]

		c = cmd.get_coords('/{}//{}/{}/CA'.format(model, chainLet[2], str(res)))[0]
		
		xyz1 = [0, 0, 0]
		xyz2 = [0, 0, 0]
		xyz3 = [0, 0, 0]

		xyz1[0] = a[0]
		xyz1[1] = a[1]
		xyz1[2] = a[2]

		xyz2[0] = b[0]
		xyz2[1] = b[1]
		xyz2[2] = b[2]

		xyz3[0] = c[0]
		xyz3[1] = c[1]
		xyz3[2] = c[2]


		def str_riffle(str_list, add):
			return(str_list[0] if len(str_list) == 1 else \
					str_list[0] + add + str_riffle(str_list[1:], add))

		def normalize(xyz):
			return xyz / np.linalg.norm(xyz)

		def get_angle(xyz1, xyz2):
			return acos(np.dot(normalize(xyz1), normalize(xyz2))) * 180 / pi

		def get_n_vec(xyz1, xyz2, xyz3):
			return normalize(np.cross(np.subtract(xyz1, xyz2), np.subtract(xyz2, xyz3)))


		cmd.reinitialize()

		#model="2bl2"
		cmd.fetch(model, async=0)

		#cmd.copy('{}_save'.format(model), model)
		
		cmd.remove('not segi A+B+C')

		cmd.color("red", '/{}/A+B+C'.format(model))

	
		nvector1=np.subtract(xyz2,xyz1)
		nvector2=np.subtract(xyz3,xyz2)
		nvector3=get_n_vec(xyz1,xyz2,xyz3)

		DotPl1=np.add(xyz1,np.multiply(nvector1,0.5))
		DotPl2=np.add(xyz2,np.multiply(nvector2,0.5))
		DotPl3=xyz2

		d1=np.dot(DotPl1,nvector1)
		d2=np.dot(DotPl2,nvector2)
		d3=np.dot(DotPl3,nvector3)

		D=np.matrix([nvector1,nvector2,nvector3])
		A=np.matrix([nvector1,nvector2,nvector3])
		A[0,0]=d1
		A[1,0]=d2
		A[2,0]=d3
		B=np.matrix([nvector1,nvector2,nvector3])
		B[0,1]=d1
		B[1,1]=d2
		B[2,1]=d3
		C=np.matrix([nvector1,nvector2,nvector3])
		C[0,2]=d1
		C[1,2]=d2
		C[2,2]=d3

		x=np.linalg.det(A)/np.linalg.det(D)
		y=np.linalg.det(B)/np.linalg.det(D)
		z=np.linalg.det(C)/np.linalg.det(D)
		point=[x,y,z]
		#print('point', point)

		cmd.pseudoatom('center', pos=[x,y,z], color='white')


		vec=normalize(np.cross(nvector1,nvector2))
		axis=[0,0,0]
		axis[0]=vec[0]
		axis[1]=vec[1]
		axis[2]=vec[2]
		#print(axis)


		ang=get_angle(nvector1,nvector2)
		n=360/ang
		nround=int(round(n))
		d2n=abs(n-nround)
		anground=360/nround
		d2ang=abs(anground-ang)

		#print('ang=',anground,'+-',d2ang)
		#print('n=',nround,'+-',d2n)

		#real_angle=anground
		#pattern = "C"
		#sel_pattern = '/{}//{}'.format(model, pattern)
		#chain_num = nround-3
		#new_models = ['{}_cp_{}'.format(model, i) for i in range(chain_num)]
		#for i in range(chain_num):
	#		cmd.create(new_models[i], sel_pattern)
	#		cmd.alter('/{}//{}'.format(new_models[i], pattern), 'chain="{}{}"'.format(pattern, i))
	#		cmd.rotate(axis, (i+1) * real_angle, '/{}//{}{}'.format(new_models[i], pattern, i), \
	#			camera=0, origin=point)
#
#		cmd.create('{}_mod'.format(model), str_riffle(['{}_cp_{}'.format(model, i) for i in range(chain_num)], ' | ') + ' | {}'.format(model))
#		for i in range(chain_num):
#			cmd.delete('{}_cp_{}'.format(model, i))


		#print(cmd.rms_cur('{}_save'.format(model), '{}_mod'.format(model)))

		
		#Mean axis and point calculation

		p1Sum =  p1Sum + point[0]
		p2Sum = p2Sum + point[1]
		p3Sum = p3Sum + point[2]


		e1MeanSum = e1MeanSum + axis[0]
		e2MeanSum = e2MeanSum + axis[1]
		e3MeanSum = e3MeanSum + axis[2]

		#print ('huinya', axis[0], e1MeanSum)

		#Mean angle and angle_error calculation
		
		
		angM = angM + ang
		sumErr = sumErr + (ang - 23.9781164782)*(ang - 23.9781164782)

	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')
	print ('')


	Mean_point[0] = p1Sum / (int(res2)-int(res1))
	Mean_point[1] = p2Sum / (int(res2)-int(res1))
	Mean_point[2] = p3Sum / (int(res2)-int(res1))

	


	Mean_axis[0] = e1MeanSum / (int(res2)-int(res1))
	Mean_axis[1] = e2MeanSum / (int(res2)-int(res1))
	Mean_axis[2] = e3MeanSum / (int(res2)-int(res1))

	angle_MEAN = angM / (int(res2)-int(res1))

	angleRMSD = sqrt(sumErr /  (int(res2)-int(res1)))

	print (angle_MEAN, '+-', angleRMSD)
	print ('vector Mean', Mean_axis)
	print ('point mean', Mean_point)







	real_angle=anground
	pattern = "C"
	sel_pattern = '/{}//{}'.format(model, pattern)
	chain_num = nround-3
	new_models = ['{}_cp_{}'.format(model, i) for i in range(chain_num)]
	for i in range(chain_num):
		cmd.create(new_models[i], sel_pattern)
		cmd.alter('/{}//{}'.format(new_models[i], pattern), 'chain="{}{}"'.format(pattern, i))
		cmd.rotate(Mean_axis, (i+1) * real_angle, '/{}//{}{}'.format(new_models[i], pattern, i), \
			camera=0, origin=Mean_point)

	cmd.create('{}_mod'.format(model), str_riffle(['{}_cp_{}'.format(model, i) for i in range(chain_num)], ' | ') + ' | {}'.format(model))
	for i in range(chain_num):
		cmd.delete('{}_cp_{}'.format(model, i))


		#alter (chain C10), chain='V'


	cmd.delete('{}'.format(model))
	cmd.fetch(model, async=0)

	
	#remove /{}_mod/A+B+C.format(model)

	#alter (segi), segi='{}'.format(chainLet[i+3])
	#print(cmd.rms_cur('{}'.format(model), '{}_mod'.format(model)))
	
			
cmd.extend('get',get_vec_point)

