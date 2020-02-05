from math import sqrt
from numpy import transpose

def sm(string):
	return (('' if string[0] != '-' else '\\') + string)

def hw5(models, base_model):
	resis_base = []
	diff = []

	cmd.iterate_state(1, '/{}////CA'.format(base_model), \
			'resis_base.append(sm(resi))', \
			space={'resis_base': resis_base, 'sm': sm})

	for m in models:
		resis = []
		cmd.iterate_state(1, '/{}////CA'.format(m), \
				'resis.append(sm(resi))', \
				space={'resis': resis, 'sm': sm})
		
		print(m, resis[1], resis[-1])
		diff.append([])

		for r in resis_base:
			ca_base = '{} and resi {} and name CA'.format(base_model, r)
			diff[-1].append(10000)
			
			for res in resis:
				ca = '{} and resi {} and name CA'.format(m, res)

				diff[-1][-1] = min(diff[-1][-1], cmd.get_distance(ca, ca_base))
	
	diff = [list(x) for x in transpose(diff)]
	dispersions = [sqrt(sum(y*y for y in x) / len(x)) for x in diff]
	
	for i in range(len(dispersions)):
		cmd.alter('/{}///{}/CA'.format(base_model, resis_base[i]), 'b="{}"'.format(dispersions[i]))

cmd.extend('hw5', hw5)