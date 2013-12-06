from math import *
from numpy import *
from myarray import *
from copy import *
from sets import Set
import sys, os
from bezier_utility import *
	   
def parse_svgfile(svg_file):
	'''
	given a svg file path, return the size of canvas,
	a list of paths.
	'''
	if not os.path.exists( svg_file ):
		raise ImportError, "SVG file not found: " + svg_file
	if not svg_file.endswith('.svg'):
		raise ImportError, svg_file + "is not SVG file."
		
	from xml.dom import minidom
	doc = minidom.parse(svg_file)
	
	svg = doc.getElementsByTagName('svg')[0]
	viewbox_string = svg.getAttribute('viewBox')
	width = float( str(svg.getAttribute('width')).replace('px', '') )
	height = float( str(svg.getAttribute('height')).replace('px', '') )
	
	view_box = str(viewbox_string).split(' ')
	view_box = [float(ele) for ele in view_box]
	assert len(view_box) == 4
	scales = [ view_box[0], view_box[1], (view_box[2] - view_box[0])/width, 
			(view_box[3] - view_box[1])/height]

	path_strings = [path.getAttribute('d') for path
				in doc.getElementsByTagName('path')]
	
 	path_strings = [str(path) for path in path_strings]
	
	doc.unlink()
	

	
	## using regex to match all the command parameters in keys
	## and match all the digits to data
	import re
	pattern = '[a-zA-Z]|-?[0-9]+\.?[0-9]*'
	matches = [re.findall(pattern, path) for path in path_strings]
	a = []
	for each in matches:
		a += each
	matches = a
	
	## test if path has any non-predefined command
	recognizable = Set(['M', 'c', 'C', 's', 'S', 'l', 'L', 'z', 'Z', 'H', 'h', 'v', 'V'])
	pattern = '[a-zA-Z]'
	keys = [re.findall(pattern, path) for path in path_strings]
	key_set = []	
	for each in keys:
		key_set += each
	key_set = Set(key_set)
	unrecognized = key_set - recognizable
	if len(unrecognized) != 0:
		raise ImportError, [x for x in unrecognized] + "cannot be recognized."
	
	## parse and samples each path	
	anchors = []
	cps = []
	samples = []
	key = None
	last_key = None
	pos = 0
	current = zeros(2)

	while(pos < len(matches)):

		if matches[pos] in recognizable: 
			last_key = deepcopy(key)
			key = matches[pos]
			pos += 1
			if (key == 'z' or key == 'Z'):
				current = asarray(deepcopy(samples[-1][:2]))
				samples[-1] += samples[-1][:2]
				last_key = key
				continue

		if (key == 'M' or key == 'm'):
			if last_key == 'M': 
				key = 'L'
				continue
			elif last_key == 'm': 
				key = 'l'
				continue
			samples.append([])
# 			if len(samples) == 8: debugger()
			cps.append([])
			anchors.append([])
			current = zeros(2)
			pts = asarray([float(x) for x in matches[pos:pos+2]])
			if key == 'm':
				cps[-1].append((pts + current).tolist())
				current += pts
			else:
				cps[-1].append(pts.tolist())
				current = pts
				
			samples[-1] += list(pts)
			anchors[-1].append(deepcopy(current))		
			pos += 2
			last_key = key
			
		elif (key == 'c' or key == 'C'):
			controls = []
			pts = asarray([float(x) for x in matches[pos:pos+6]]).reshape(-1,2)
			if key == 'c':
				controls = [deepcopy(current).tolist()]+(pts + current).tolist()
				cps[-1] += (pts + current).tolist()
				current += pts[-1]
			else: 
				controls = [deepcopy(current).tolist()]+pts.tolist()
				cps[-1] += pts.tolist()
				current = pts[-1]
			controls = asarray(controls).reshape(4,2)
			samples[-1] += sample_cubic_bezier_curve(controls)[0][2:]
			
			anchors[-1].append(deepcopy(current))		
			pos += 6
			last_key = key
			
		elif (key == 's' or key == 'S'):
			if last_key in Set(['S', 's', 'C', 'c']):
				ctrl1 = current*2 - asarray(cps[-1][-2])
			else:
				ctrl1 = current
			controls = []
			pts = asarray([float(x) for x in matches[pos:pos+4]]).reshape(-1,2)
			if key == 's':	
				controls = [deepcopy(current).tolist()]+[ctrl1.tolist()]+(
							pts + current).tolist()
				cps[-1] += [ctrl1.tolist()]+(pts + current).tolist()
				current += pts[-1]
			else:
				controls = [deepcopy(current).tolist()]+[ctrl1.tolist()]+pts.tolist()	
				cps[-1] = [ctrl1.tolist()]+pts.tolist()
				current = pts[-1]
							
			controls = asarray(controls).reshape(4,2)
			samples[-1] += sample_cubic_bezier_curve(controls)[0][2:]
			
			anchors[-1].append(deepcopy(current))
			pos += 4
			last_key = key

		elif (key == 'l' or key == 'L'):
			begin = current[:]
			pts = asarray([float(x) for x in matches[pos:pos+2]])
			if key == 'l':
				cps[-1] += (pts + current).tolist()
				current += pts
			else: 
				cps[-1] += pts.tolist()
				current = pts
			end = current[:]
			samples[-1] += sample_straight_line(begin, end)[0][2:]
			anchors[-1].append(deepcopy(current))
			pos += 2
			last_key = key
			
		elif (key == 'h' or key == 'H'):
			begin = current[:]
			pts = float(matches[pos])
			if key == 'h':
				cps[-1] += [current[0]+pts, current[1]]
				current[0] += pts
			else: 
				cps[-1] += [pts, current[1]]
				current[0] = pts
			end = current[:]
			samples[-1] += sample_straight_line(begin, end)[0][2:]
			anchors[-1].append(deepcopy(current))
			pos += 1
			last_key = key
			
		elif (key == 'v' or key == 'V'):
			begin = current[:]
			pts = float(matches[pos])
			if key == 'v':
				cps[-1] += [current[0], current[1]+pts]
				current[1] += pts
			else: 
				cps[-1] += [current[0], pts]
				current[1] = pts
			
			end = current[:]
			samples[-1] += sample_straight_line(begin, end)[0][2:]
			anchors[-1].append(deepcopy(current))
			pos += 1
			last_key = key

	for curve in samples:
		curve[::2] = ((asarray(curve[::2])-scales[0])/scales[2]).tolist()
		curve[1::2] = ((asarray(curve[1::2])-scales[0])/scales[2]).tolist()

	for i, each in enumerate(samples):
		if len(each) < 4:
			print i, each

	return [width, height], samples
	
if __name__ == '__main__': main()
	