from pymol import cmd
import numpy as np

def sm(string):
	return (('' if string[0] != '-' else '\\') + string)

def hw3(model):
	
	cmd.remove(model + ' and resn hoh') 
	cmd.remove(model + ' and not alt a+""')
	cmd.remove(model + ' and not chain A')


	resis = []
	cmd.iterate('/{}//A//CA'.format(model), 'resis.append(resi)', space={'resis':resis})
	
	resis_int = map(int, resis)
	resis = map(sm, resis)
	
		
	
	#Pro

	resP = []
	cmd.iterate('/{}//A//CA'.format(model), 'resP.append(resn == "PRO")', space={'resP':resP})
	
	
	angs_phi = []
	angs_psi = []
	for i in range(1, len(resis) - 1):
		if (resis_int[i-1] + 1 == resis_int[i]) and (resis_int[i] + 1 == resis_int[i+1]) and resP[i]:
			angs_phi.append(cmd.get_dihedral('/{}//A/{}/C'.format(model, resis[i-1]), \
					'/{}//A/{}/N'.format(model, resis[i]), \
					'/{}//A/{}/CA'.format(model, resis[i]), \
					'/{}//A/{}/C'.format(model, resis[i]))),
			angs_psi.append(cmd.get_dihedral('/{}//A/{}/N'.format(model, resis[i]), \
					'/{}//A/{}/CA'.format(model, resis[i]), \
					'/{}//A/{}/C'.format(model, resis[i]), \
					'/{}//A/{}/N'.format(model, resis[i+1])))
	fname = 'C:\\Users\\Nick\\Desktop\\Bank\\PhystechStuding\\Buslaev\\project3\\hw3_{}_Pro.dat'.format(model)
	open(fname, 'w').close()
	np.savetxt(fname, np.transpose([angs_phi, angs_psi]))

cmd.extend('hw3', hw3)