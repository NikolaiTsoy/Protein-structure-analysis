from pymol import cmd
import numpy as np
from math import acos, pi
from statistics import pstdev, mean


mod_string = 'mod'


def str_riffle(str_list, add=' '):
	return '' if len(str_list) == 0 else \
			(str_list[0] if len(str_list) == 1 else \
			str_list[0] + add + str_riffle(str_list[1:], add))

def normalize(xyz):
	return xyz / np.linalg.norm(xyz)

def get_n_vec(xyz1, xyz2, xyz3):
	return normalize(np.cross(np.subtract(xyz1, xyz2), np.subtract(xyz2, xyz3)))

def get_angle(xyz1, xyz2):
	return acos(np.dot(xyz1, xyz2)) * 180 / pi

def new_segi(segi, all_segis):
	'''
	Make new segi name from old

	If segi model names [A,B,C,D,E], used names [A,B,C], then created will be:
		[D,E,Fmod,Gmod,Hmod,Imod,Jmod]
	but we need:
		[D,E,Amod,Bmod,Cmod,Dmod,...]
	'''
	return chr(ord(all_segis[0]) + (ord(segi[0])-ord(all_segis[0])) % len(all_segis)) + \
			(segi[1:] if len(segi) > 1 else '')


def find_lines(xyzs):
	z = 0
	ns = []
	Ds = []
	for i in range(len(xyzs) - 1):
		# (x, n) + D = 0
		# A*x + B*y + C*z + D = 0
		ns.append(np.subtract(xyzs[i], xyzs[i+1]))
		point = [(xyzs[i][j] + xyzs[i+1][j]) / 2 for j in range(3)]
		Ds.append(-np.dot(point, ns[-1]))

	line_points = []
	line_ns = []
	for i in range(len(xyzs) - 2):
		# plane intersection point:
		# 	(-B2*(D1+C1*z) + B1*(D2+C2*z)) / (A1*B2-A2*B1),
		# 	( A2*(D1+C1*z) - A1*(D2+C2*z)) / (A1*B2-A2*B1),
		# 	z
		# line:
		# point + lambda * normalize([n1 x n2])
		det = ns[i][0]*ns[i+1][1] - ns[i][1]*ns[i+1][0]
		
		line_points.append([\
				(-ns[i+1][1]*(Ds[i]+ns[i][2]*z) + ns[i][1]*(Ds[i+1]+ns[i+1][2]*z)) / det, \
				( ns[i+1][0]*(Ds[i]+ns[i][2]*z) - ns[i][0]*(Ds[i+1]+ns[i+1][2]*z)) / det, \
				z \
			])

		line_ns.append(normalize(np.cross(ns[i], ns[i+1])))

	return line_points, line_ns


def find_mean_line(xyzs_resis):
	line_points, line_ns = np.transpose([find_lines(xyzs) for xyzs in xyzs_resis], \
			axes=[1,0,2,3])
	line_points =[e1 for e2 in line_points for e1 in e2]
	line_ns = [e1 for e2 in line_ns for e1 in e2]

	n = np.mean(line_ns, axis=0).tolist()
	point = np.mean(line_points, axis=0).tolist()

	n_error = [pstdev(line_n[i] for line_n in line_ns) for i in range(3)]
	point_error = [pstdev(line_point[i] for line_point in line_points) for i in range(3)]

	return n, n_error, point, point_error


def make_segis(model, resi1, resi2, resi3, absent_segis_names):
	segis = []
	for resi in (resi1, resi2, resi3):
		segis.append([])
		cmd.iterate('/{}///{}/CA'.format(model, resi), \
				'segis[-1].append(segi)', \
				space={'segis': segis})
	segis = [set(x) for x in segis]
	segis = sorted(segis[0] & segis[1] & segis[2])

	# xyzs_segis format:
	# [
	# 	{segi1}[
	# 		{resi1}[x,y,z],
	# 		{resi2}[x,y,z],
	# 		{resi3}[x,y,z]
	# 	],
	# 	{segi2}[...],
	# 	{segi3}[...],
	# 	... 	
	# ]
	xyzs_segis = \
		cmd.get_coords('/{}/{}//{}+{}+{}/CA'.\
			format(model, str_riffle(segis, '+'), resi1, resi2, resi3))
	xyzs_segis = np.split(xyzs_segis, len(segis))
	
	# normal vectors for each segi
	n_vecs = [get_n_vec(*xyzs_segi) for xyzs_segi in xyzs_segis]
	# angles between each adjacent pair
	angles = [get_angle(n_vecs[i], n_vecs[i + 1]) for i in range(len(n_vecs) - 1)]

	# rotation angle between two adjacent aminoacids with error
	angle_value, angle_error = mean(angles), pstdev(angles)
	segi_count = int(round(360/angle_value))
	real_angle = 360 // segi_count

	print(angle_value, angle_error)
	print(real_angle)

	# main axis parameters
	axis, axis_error, point, point_error = \
			find_mean_line(np.transpose(xyzs_segis, axes=[1,0,2]))

	# create missing segis (by creating objects) by copying last segi
	# if segi model names [A,B,C,D,E], used names [A,B,C], then created will be:
	# 	[D,E,Fmod1,Gmod1,Hmod1,Imod1,Jmod1,Kmod2,Lmod2,...]
	#
	# rotate objects using parameters found previously
	pattern = segis[-1]
	sel_pattern = '/{}/{}'.format(model, pattern)
	segi_num = int(round(360 / angle_value - len(segis)))
	new_models = ['{}_cp_{}'.format(model, i) for i in range(segi_num)]
	for i in range(segi_num):
		segi_name = absent_segis_names[(i+len(segis)) % len(absent_segis_names)]
		new_segi_count = (i+len(segis)) // len(absent_segis_names)
		mod_string_i = '{}{}'.format(mod_string, new_segi_count) if new_segi_count != 0 else ''

		cmd.create(new_models[i], sel_pattern)
		cmd.alter('/{}/{}'.format(new_models[i], pattern), \
				'segi="{}{}"'.format(segi_name, mod_string_i))
		cmd.rotate(axis, (i+1) * real_angle, \
				'/{}/{}{}'.format(new_models[i], segi_name, mod_string_i), \
				camera=0, origin=point)

	### CREATE FINAL OBJECT by merging accessory models ###
	new = '{}_calc'.format(model)
	cmd.create(new, \
			(str_riffle(new_models, ' | ') + ' | ' if len(new_models) != 0 else '') + \
			'/{}/{}'.format(model, str_riffle(segis, '+')))
	
	# delete accessory models
	for new_model in new_models:
		cmd.delete(new_model)

	return new, segi_count

def merge(fin, models):
	cmd.create('{}_sym'.format(fin), str_riffle(models, ' | '))
	for model in models:
		cmd.delete(model)

def prepare(model, resi1, resi2, resi3, segis_for_save, segis_for_del):
	resi1,resi2,resi3 = sorted((resi1,resi2,resi3))
	
	cmd.fetch(model, async=0)
	cmd.remove('not alt a+""')

	# delete additional protein segis
	if segis_for_del != '':
		cmd.remove('segi {}'.format(segis_for_del.replace(' ', '+')))

	# all_segis contains segis which are included to that protein part (symmetric) we need
	all_segis = []
	for resi in (resi1, resi2, resi3):
		all_segis.append([])
		cmd.iterate('/{}///{}/CA'.format(model, resi), \
				'all_segis[-1].append(segi)', \
				space={'all_segis': all_segis})
	all_segis = [set(x) for x in all_segis]
	all_segis = all_segis[0] | all_segis[1] | all_segis[2]

	cmd.remove('not segi {}'.format(str_riffle(list(all_segis), '+')))

	if type(segis_for_save) == list:
		segis_for_save_fact = set(segis_for_save)
		absent_segis_names = sorted(all_segis - segis_for_save_fact)
	else:
		segis_for_save_fact = absent_segis_names = sorted(all_segis)

	model_save = '{}_save'.format(model)
	cmd.copy(model_save, model)
	cmd.remove('{} and not segi {}'.format(model, str_riffle(list(segis_for_save_fact), '+')))

	return sorted(all_segis)

def project(model, resi1_str, resi2_str, resi3_str, \
		symQ, syms_for_use, syms_for_del, sym_d=3, \
		segis_for_save='all', segis_for_del=''):
	all_segis = prepare(model, resi1_str, resi2_str, resi3_str, segis_for_save, segis_for_del)
	
	new_model, segi_count = make_segis(model, resi1, resi2, resi3, all_segis)
	cmd.alter(new_model, \
			'segi = new_segi(segi, sorted(all_segis))', \
			space={'new_segi': new_segi, 'all_segis': all_segis, 'sorted': sorted})

	if symQ:
		cmd.symexp('{}_sym'.format(model), \
				'{}_save'.format(model), \
				'({}_save)'.format(model), \
				sym_d)

		for i in range(len(syms_for_use)):
			cmd.alter('{}_sym{}'.format(model, syms_for_use[i]), 'segi=segi+"mod{}"'.format(i + 1))

		merge(model, ['{}_save'.format(model)] + ['{}_sym{}'.format(model, s) for s in syms_for_use])
	
	else: merge(model, ['{}_save'.format(model)])
	
	for sym in syms_for_del:
		cmd.delete('{}_sym{}'.format(model, sym))
	
	cmd.alter('all', 'chain="A"')

	cmd.save('{}_calc.cif'.format(model), '{}_calc'.format(model))
	print(cmd.rms_cur('{}_calc'.format(model), '{}_sym'.format(model)))

cmd.extend('project', project)