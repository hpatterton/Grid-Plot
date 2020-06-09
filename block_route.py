import graphics as gfx
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# import os


class Graphics_Svg():

	def __init__(self, filename):
		self.outfilename = filename
		ini_file = open('svg.ini','r')
		self.svg_text = ini_file.read()
		self.object_id = 0

	def draw_rect_svg(self, x1, y1, width, height, fill_colour, stroke_colour, stroke_width):
		self.svg_text += '\t\t<rect\n'
		self.svg_text += '\t\tstyle = \"fill:#'+fill_colour+';fill-opacity:1;stroke:#'+stroke_colour+';stroke-width:'+str(stroke_width)+';stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0\"\n'
		self.svg_text += '\t\tid = \"'+str(self.object_id)+'\"\n'
		self.svg_text += '\t\twidth = \"'+str(width)+'\"\n'
		self.svg_text += '\t\theight = \"'+str(height)+'\"\n'
		self.svg_text += '\t\tx = \"'+str(x1)+'\"\n'
		self.svg_text += '\t\ty = \"'+str(y1)+'\"/>\n'
		self.object_id += 1

	def draw_circle_svg(self, x1, y1, radius, fill_colour, stroke_colour, stroke_width):
		self.svg_text += '\t\t<circle\n'
		self.svg_text += '\t\tstyle = \"fill:#'+fill_colour+';fill-opacity:1;stroke:#'+stroke_colour+';stroke-width:'+str(stroke_width)+';stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0\"\n'
		self.svg_text += '\t\tid = \"'+str(self.object_id)+'\"\n'
		self.svg_text += '\t\tcx = \"'+str(x1)+'\"\n'
		self.svg_text += '\t\tcy = \"'+str(y1)+'\"\n'
		self.svg_text += '\t\tr = \"' + str(radius) + '\"/>\n'
		self.object_id += 1

	def draw_line_svg(self, x1, y1, x2, y2, colour, stroke_width):
		self.svg_text += '\t\t<path\n'
		self.svg_text += '\t\tstyle = \"fill:none;stroke:#'+str(colour)+';stroke-width:'+str(stroke_width)+'px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n'
		self.svg_text += '\t\td = \"M '+str(x1)+','+str(y1)+' L '+str(x2)+','+str(y2)+'\"\n'
		self.svg_text += '/>\n'

	def draw_labels_svg(self,label_string, x1, y1, font, size, colour, horizontal_alignment):
		self.svg_text += '\t\t<text\n'
		self.svg_text += '\t\tx = \"'+str(x1)+'\"\n'
		self.svg_text += '\t\ty = \"'+str(y1)+'\"\n'
		self.svg_text += '\t\tstyle = \"font-family: '+font+';\n'
		self.svg_text += '\t\ttext-anchor: '+horizontal_alignment+';\n'
		self.svg_text += '\t\tfont-size: '+str(size)+'px;\n'
		self.svg_text += '\t\tstroke:none;stroke-width:0.26458332\n'
		self.svg_text += '\t\tfill:#ffff00;\n'
		self.svg_text += '\t\tfill-opacity:1;\"'
		self.svg_text += '\t\t>'+label_string+'</text>\n'


	def write_file(self):
		outfile = open(self.outfilename,'w')
		self.svg_text += '\t\t</g>\n'
		self.svg_text += '\t\t</svg>\n'
		outfile.write(self.svg_text)
		outfile.close()

class Block_Route:

	def __init__(self, list_of_categories=0):
		#self.write_to_file = write_to_file
		self.longest_path = 0
		self.data_array = []
		self.number_of_entries = 0
		# the order of the categories must match the index used in the data file
		if (list_of_categories == 0):
			self.category = ['Inversion', 'Translocation (balanced)', 'Translocation (unbalanced)', 'Fusion', 'Fission',
							 'Transposition', 'block interchange']
		else:
			self.category = list_of_categories
		self.number_of_categories = len(self.category)
		self.summed_array = []

	def ReadFile(self,filepath):
		handle = open(filepath,'r')
		for line in handle:
			temp = line.rstrip('\n').split(',')
			self.data_array.append(temp)
			if len(temp) > self.longest_path:
				self.longest_path = len(temp)
		self.number_of_entries = len(self.data_array)

	def GetLongestPath(self):
		return self.longest_path

	def GetNumberOfCategories(self):
		return self.number_of_categories

	def GenerateSummedTransitionsMatrixOneRun(self):
		number_of_array_columns = self.GetLongestPath()
		number_of_array_rows = len(self.data_array)
		summed_array_row = [0,]*number_of_array_columns
		for level in range(number_of_array_columns):
			summed_level = []
			for row in range(number_of_array_columns):
				summed_level.append(summed_array_row.copy())
			self.summed_array.append(summed_level.copy())

		for column in range(number_of_array_columns):
			for row in range(number_of_array_rows):
				if column < len(self.data_array[row]) - 1:
					first_entry = int(self.data_array[row][column])
					second_entry = int(self.data_array[row][column+1])
					self.summed_array[column][first_entry][second_entry] += 1
		#self.PrintSummedMatrix(self.summed_array)
		return self.summed_array

	def NormalizeSummedMatrix(self, normalize_global=False, maximum_value=5):

		"""This function finds the maximum for all transitions in a run
		or the maximum for a specific transition (per level), eg. column 0 to 1, or
		1 to 2, etc.  All values (global or per column pair transition)
		 are normalized to the respective maximum (or maxima), and adjusted to a
		 maximum value specified by the user."""

		normalized_data_matrix = []
		if normalize_global==False:
			normalized_array = []
			for level in range(len(self.summed_array)):
				maximum = 0
				for row in range(len(self.summed_array[level])):
					for column in range(len(self.summed_array[level][row])):
						if self.summed_array[level][row][column] > maximum:
							maximum = self.summed_array[level][row][column]
				rows = []
				for row in range(len(self.summed_array[level])):
					local_row = [0, ] * len(self.summed_array[level][row])
					for column in range(len(self.summed_array[level][row])):
						if maximum > 0:
							local_row[column] = maximum_value * self.summed_array[level][row][column] / maximum
						else:
							local_row[column] = self.summed_array[level][row][column]
					rows.append(local_row)
				normalized_array.append(rows.copy())
		else:
			maximum = 0
			for level in range(len(self.summed_array)):
				for row in range(len(self.summed_array[level])):
					for column in range(len(self.summed_array[level][row])):
						if self.summed_array[level][row][column] > maximum:
							maximum = self.summed_array[level][row][column]
			normalized_array = []
			h_levels=len(self.summed_array)
			for level in range(len(self.summed_array)):
				rows = []
				h_rows=len(self.summed_array[level])
				for row in range(len(self.summed_array[level])):
					h_column= len(self.summed_array[level][row])
					local_row = [0, ] * len(self.summed_array[level][row])
					for column in range(len(self.summed_array[level][row])):
						if maximum > 0:
							local_row[column] = maximum_value * self.summed_array[level][row][column] / maximum
						else:
							local_row[column] = self.summed_array[level][row][column]
					rows.append(local_row)
				normalized_array.append(rows.copy())
		return normalized_array

	def PrintSummedMatrix(self, matrix):
		for level in range(len(matrix)):
			print('Level =', level)
			for row in range(len(matrix[level])):
				print(matrix[level][row])


class Block_Plot():

	def __init__(self, datafile, svg_outfile, grid_hspace=30, grid_vspace=20, max_linewidth=5, normalize_global=True, categories=0):
		if categories==0:
			self.category = ['Inversion', 'Translocation (balanced)', 'Translocation (unbalanced)', 'Fusion', 'Fission',
						 	 'Transposition', 'Block interchange']
		else:
			self.category = categories

		self.start_columns = 70
		self.start_rows = 20

		# block rectangle
		self.top_left_x=0
		self.top_left_y=0
		self.bottom_right_x=5
		self.bottom_right_y=5

		# text settings
		self.text_pointsize = 10
		self.text_font_family ="arial"

		self.scale = 0.5
		self.inter_column_space = grid_hspace
		self.inter_row_space = grid_vspace
		self.graphics_window = gfx.GraphWin("Block Route", 2970, 2100)
		self.graphics_window.setBackground("white")

		# read and normalize data
		data = Block_Route()
		data.ReadFile(datafile)
		print(data.data_array)
		self.longest_path = data.longest_path
		self.number_of_rows = data.number_of_categories
		self.number_of_columns = data.number_of_entries
		data.GenerateSummedTransitionsMatrixOneRun()
		self.normalised_array = data.NormalizeSummedMatrix(normalize_global, max_linewidth)
		data.PrintSummedMatrix(self.normalised_array)
		self.svg_file = Graphics_Svg(svg_outfile)
		self.draw_edges()
		self.draw_vertexes()


	def draw_vertexes(self):
		for row in range(0,self.number_of_rows):
			for column in range(0,self.longest_path):
				previous_x = self.scale*(self.start_columns+column*(self.inter_column_space+self.bottom_right_x/2))
				previous_y = self.scale*(self.start_rows+row*(self.inter_row_space+self.bottom_right_y/2))
				new_point_x = previous_x+self.scale*(self.inter_column_space+self.bottom_right_x/2)
				new_point_y = previous_y+self.scale*(self.inter_row_space+self.bottom_right_y/2)
				#bottom_right_point_x = previous_x+self.scale*(self.inter_column_space+self.bottom_right_x)
				#bottom_right_point_y = previous_y+self.scale*(self.inter_row_space+self.bottom_right_y)
				#r = gfx.Rectangle(gfx.Point(new_point_x, new_point_y), gfx.Point(new_point_x, bottom_right_point_y))
				self.svg_file.draw_circle_svg(new_point_x, new_point_y, self.scale*self.bottom_right_y/2, 'ffffff', '000000', 0.5)
				#self.svg_file.draw_rect_svg(top_left_point_x, top_left_point_y, bottom_right_point_x-top_left_point_x, bottom_right_point_y-top_left_point_y, 'ffffff', '000000', 0.5)

				#r.setFill("white")
				#r.draw(self.graphics_window)

		for row in range(0, self.number_of_rows):
			#x_value = self.scale * self.start_columns / 1.5
			y_value = self.scale * self.start_rows + self.scale * row * ( self.inter_row_space+ self.bottom_right_y/2) + self.scale * (
						self.inter_row_space + self.bottom_right_y / 2)
			#t = gfx.Text(gfx.Point(x_value, y_value), self.category[row])
			#t.setSize(self.text_pointsize)
			#t.setFace(self.text_font_family)
			#t.draw(self.graphics_window)
			x_value = self.scale*(self.start_columns+self.inter_column_space)-8
			self.svg_file.draw_labels_svg(self.category[row], x_value, y_value+1, 'Arial', 6, '000000', 'end')

		self.svg_file.write_file()


	def draw_edges(self):
		"""This function receives a stack of 2D arrays (list of list of list).
		Each individual array is connections between blocks n column i and i+1
		and each stack is i from 0 to N-1 for N columns in block plot"""
		# categories = len(self.category)
		# weight_array = []
		# columns = [1,1,2,3,4,5,6]
		# for i in range(categories):
		# 	weight_array.append(columns)
		for column_left in range(0,self.longest_path-1):
			for row_left in range(0, self.number_of_rows ):
				for row_right in range(0, self.number_of_rows ):
					start_x = self.scale * (self.start_columns + (column_left+1) * (self.inter_column_space + self.bottom_right_x/2))
					start_y = self.scale * (self.start_rows + (row_left+1) * (self.inter_row_space + self.bottom_right_y/2))
					end_x = start_x+self.scale *  (self.inter_column_space + self.top_left_x+self.bottom_right_x/2)
					end_y = self.scale * (self.start_rows + (row_right+1) * (self.inter_row_space +  self.bottom_right_y/2))
					#line = gfx.Line(gfx.Point(start_x, start_y), gfx.Point(end_x, end_y))
					#line.setWidth(1*weight_array[row_left][row_right])
					#line.draw(self.graphics_window)
					#print(column_left, row_left, row_right)
					self.svg_file.draw_line_svg(start_x, start_y, end_x, end_y, '000000', 0.5*self.normalised_array[column_left][row_left][row_right])







