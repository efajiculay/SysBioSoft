from pylint.lint import Run
import os

summary = {}
for dirpath, dirs, files in os.walk('BioSANS2020'): 
	for filename in files:
		fname = os.path.join(dirpath,filename)
		if fname.endswith('.py'):
			results = Run([fname], do_exit=False)
			try:
				summary[filename] = results.linter.stats['global_note']
			except:
				pass

print()
for x in summary:			
	print(x,summary[x])

   