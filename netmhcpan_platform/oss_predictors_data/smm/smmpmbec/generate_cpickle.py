import os
import re
import csv
import cPickle


def run():
	dir_name = 'new_trained_data'
	output = 'output' 
	path_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), dir_name)
	for file in os.listdir(dir_name):
		m = re.match(r'(?!.*model)(.*)\.txt', file)
		if m:
			path_file = os.path.join(path_dir, file)
			data = csv.reader(open(path_file), delimiter='\t')
			length = data.next()[1]
			offset = ''
			mat = {}
			for row in data:
				if row:
					if row[0] == 'Intercept':
						intercept = row[1]
					else: 
						scores = [float(i) for i in row[1:]]
						mat[row[0]] = tuple(scores)
						
			filename = m.group(1).replace(':','')
			filename = os.path.join(output, filename)+'.cpickle'
			try:
    				os.remove(filename)
			except OSError:
    				pass
			fout = open(filename,"wb")
			cPickle.dump(int(length), fout)
			cPickle.dump(mat, fout)
			cPickle.dump(float(intercept), fout)
			fout.close()
	print "All done!"

#if __name__ == '_main __':
run()

