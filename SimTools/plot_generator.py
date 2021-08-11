import matplotlib.pyplot as plt
import numpy as np
import random
import time
import datetime
import math

class PlotGenerator:
	"""
	Class for generating plots using matplotlib.
	"""
	def __init__(self, use_random=False, reproducability_seed=None, number_of_colors=16, color_difference_threshold=0.20):
		"""
		Initialize self.
		Sets all settables to their default values.
		WARNING: setting the number_of_colors and/or color_difference_threshold to large will cause long if not infinite wait times for color generation.
		"""
		self.__x_scale = None #matplotlib defaults to linear if not specified
		self.__x_base = 10
		self.__y_scale = None #matplotlib defaults to linear if not specified
		self.__y_base = 10
		self.__x_tick_label_rotation = 0
		self.__x_tick_label_horizontal_alignment = 'center'
		self.__y_tick_label_rotation = 0
		self.__y_tick_label_horizontal_alignment = 'center'
		self.__title = ""
		self.__x_label = ""
		self.__y_label = ""
		self.__x_tick_labels = None
		self.__y_tick_labels = None
		self.__x_tick_label_font_size = 12
		self.__figure_width = 6.4
		self.__figure_height = 4.8
		self.__COLORS = ['r', 'b', 'g', 'y', 'm', 'c'] #Use matplotlib's "Base Colors"
		if use_random: #Use random but non-similar colors 
			seed = reproducability_seed if reproducability_seed != None else self.generate_seed()
			random.seed(seed)
			print('Seed:', seed)
			self.__COLORS = self.generate_colors(number_of_colors, color_difference_threshold)

	def set_figure_options(self, width=6.4, height=4.8):
		self.__figure_width = width
		self.__figure_height = height

	def set_text(self, title, x_label, y_label, x_tick_labels, y_tick_labels):
		"""
		Sets the titles and labels of to-be-generated plots.
		"""
		if x_tick_labels != None:
			assert isinstance(x_tick_labels, list), "x_tick_labels must be a list"
		if y_tick_labels != None:
			assert isinstance(y_tick_labels, list), "y_tick_labels must be a list"
		self.__title = title
		self.__x_label = x_label
		self.__y_label = y_label
		self.__x_tick_labels = x_tick_labels
		self.__y_tick_labels = y_tick_labels

	def set_text_options(self, x_tick_label_rotation, x_tick_label_horizontal_alignment, y_tick_label_rotation, y_tick_label_horizontal_alignment, x_tick_label_font_size):
		"""
		Sets the title and label text options of to-be-generated plots.
		"""
		self.__x_tick_label_rotation = x_tick_label_rotation
		self.__x_tick_label_horizontal_alignment = x_tick_label_horizontal_alignment
		self.__y_tick_label_rotation = y_tick_label_rotation
		self.__y_tick_label_horizontal_alignment = y_tick_label_horizontal_alignment
		self.__x_tick_label_font_size = x_tick_label_font_size

	def set_axis_options(self, x_scale, x_base, y_scale, y_base):
		"""
		Set the axis options of to-be-generated plots.
		Changing the bases does nothing if the scale is None or 'linear'.
		Not all plots can change scales.
		"""
		if x_scale != 'linear':
			self.__x_scale = None
		else:
			self.__x_scale = x_scale
		self.__x_base = x_base
		if y_scale == 'linear':
			self.__y_scale = None
		else:
			self.__y_scale = y_scale
		self.__y_base = y_base

	def bar_plot(self, data_sets, data_set_names):
		"""
		Generates a bar plot of one or many data sets.
		Cannot change x axis scale.
		"""
		assert isinstance(data_sets, (list, np.ndarray)), 'data_sets must be a list or np.ndarray'
		for data_set in data_sets:
			assert isinstance(data_set, np.ndarray), 'each data set in data_sets must be a np.ndarray'
		NUM_SETS = len(data_sets)
		assert isinstance(data_set_names, list), 'data_set_names must be a list'
		assert NUM_SETS == len(data_set_names), 'data_sets length must equal data_set_names length'
		assert NUM_SETS <= len(self.__COLORS), 'data_sets length must be less than or equal to the number of colors'
		for i in range(1, NUM_SETS):
			assert len(data_sets[i]) == len(data_sets[i - 1]), 'all data sets must be of equal length'
		NUM_BARS = len(data_sets[0])
		for ds in data_sets:
			assert len(ds) == len(self.__x_tick_labels), 'all data sets must be of length equal to the number of x tick labels'

		plt.figure(figsize=(self.__figure_width, self.__figure_height))

		WIDTH = 1 / (NUM_SETS + 1)
		HALF_WIDTH = WIDTH / 2
		SHIFT_TO_CENTER = (NUM_SETS * HALF_WIDTH)
		x = np.arange(NUM_BARS)
		for i in range(0, NUM_SETS):
			y = data_sets[i]
			i_h = NUM_SETS - (i + 1) % NUM_SETS
			OFFSET = HALF_WIDTH * (NUM_SETS - 2 * ((i + 1) % NUM_SETS) - 1)
			plt.bar(x + OFFSET, y, color=self.select_color(i), width=WIDTH)

		if self.__y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.__y_scale, basey=self.__y_base)

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)
		if self.__x_tick_labels != None:
			x = np.arange(NUM_BARS)
			plt.xticks(x, labels=self.__x_tick_labels, rotation=self.__x_tick_label_rotation, ha=self.__x_tick_label_horizontal_alignment, fontsize=self.__x_tick_label_font_size)

		plt.legend(data_set_names, loc='upper left')

		plt.show()

	def box_plot(self, data_sets, data_set_names, showfliers):
		"""
		Generates a box plot of one or many data sets.
		Cannot change x axis scale.
		"""
		assert isinstance(data_sets, (list, np.ndarray)), 'data_sets must be a list or np.ndarray'
		for data_set in data_sets:
			assert isinstance(data_set, np.ndarray), 'each data set in data_sets must be a np.ndarray'
		NUM_SETS = len(data_sets)
		assert isinstance(data_set_names, list), 'data_set_names must be a list'
		assert NUM_SETS == len(data_set_names), 'data_sets length must equal data_set_names length'
		assert NUM_SETS <= len(self.__COLORS), 'data_sets length must be less than or equal to the number of colors'
		for i in range(1, NUM_SETS):
			assert len(data_sets[i]) == len(data_sets[i - 1]), 'all data sets must be of equal length'
		for ds in data_sets:
			assert len(ds) == len(self.__x_tick_labels), 'all data sets must be of length equal to the number of x tick labels'

		plt.figure(figsize=(self.__figure_width, self.__figure_height))

		boxes = []
		for i in range(0, NUM_SETS):
			boxes.append(self.__gen_box_plot_object(data_sets[i], self.select_color(i), i, NUM_SETS, showfliers))

		if self.__y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.__y_scale, basey=self.__y_base)

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)
		if self.__x_tick_labels != None:
			new_x_tick_labels = self.combine_data_set_names_with_x_tick_labels(data_set_names)
			x_ticks = np.arange(0, NUM_SETS * len(data_sets[0]), 1)
			plt.xticks(x_ticks, labels=new_x_tick_labels, rotation=self.__x_tick_label_rotation, ha=self.__x_tick_label_horizontal_alignment, fontsize=self.__x_tick_label_font_size)

		for_legend = []
		for box in boxes:
			for_legend.append(box['boxes'][0])
		plt.legend(for_legend, data_set_names, loc='upper left')

		plt.show()

	def __gen_box_plot_object(self, data, color, plot_num, num_plots, showfliers):
		"""
		Creates a box plot object.
		Used in box_plot function.
		"""
		positions = np.arange(plot_num, len(data) * num_plots, num_plots)

		box_plot = plt.boxplot(data, patch_artist=True, positions=positions, showfliers=showfliers)
		for box in box_plot['boxes']:
			box.set(color=color, linewidth=1)
			box.set(facecolor='white')
		for flier in box_plot['fliers']:
			flier.set(markersize=0.5, alpha=0.5)
		return box_plot

	def histogram(self, data_sets, bins, data_set_names):
		"""
		Generates a histogram.
		Cannot change x axis scale.
		"""
		for data in data_sets:
			assert isinstance(data, np.ndarray), 'data must be a np.ndarray'
		
		plt.figure(figsize=(self.__figure_width, self.__figure_height))

		for i in range(len(data_sets)):
			data = data_sets[i]
			plt.hist(data, bins, color=self.select_color(i))

		if self.__y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.__y_scale, basey=self.__y_base)

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)
		plt.legend(data_set_names, loc='upper left')

		plt.show()

	def scatter_plot(self, data_x, data_y, trendline=False):
		"""
		Generates a scatter plot.
		Cannot change x or y axis scale.
		"""
		assert isinstance(data_x, np.ndarray), 'data_x must be a np.ndarray'
		assert isinstance(data_y, np.ndarray), 'data_y must be a np.ndarray'

		plt.figure(figsize=(self.__figure_width, self.__figure_height))

		plt.scatter(data_x, data_y)

		if trendline:
			z = np.polyfit(data_x, data_y, 1)
			p = np.poly1d(z)
			plt.plot(data_x, p(data_x), 'r--')

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)

		plt.show()

	def heatmap(self, data_sets):
		"""
		Generates a heatmap.
		Cannot change x or y axis scale.
		"""
		assert isinstance(data_sets, (list, np.ndarray)), 'data_sets must be a list or np.ndarray'
		for data_set in data_sets:
			assert isinstance(data_set, np.ndarray), 'each data set must be a np.ndarray'
		NUM_SETS = len(data_sets)
		if NUM_SETS > 1:
			for i in range(1, NUM_SETS):
				assert len(data_sets[i]) == len(data_sets[i - 1]), 'all data sets must be of equal length'

		plt.figure(figsize=(self.__figure_width, self.__figure_height))

		plt.imshow(data_sets, cmap='hot', interpolation='nearest')

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)

		if self.__x_tick_labels != None:
			assert len(self.__x_tick_labels) == len(data_sets[0])
			plt.xticks(np.arange(len(data_sets[0])), self.__x_tick_labels, rotation=self.__x_tick_label_rotation, ha=self.__x_tick_label_horizontal_alignment, fontsize=self.__x_tick_label_font_size)
		if self.__y_tick_labels != None:
			assert len(self.__y_tick_labels) == NUM_SETS
			plt.yticks(np.arange(len(self.__y_tick_labels)), self.__y_tick_labels, rotation=self.__y_tick_label_rotation, ha=self.__y_tick_label_horizontal_alignment)

		plt.show()

	def line_plot(self, x, data_sets, trendline=False):
		"""
		Generates a line plot.
		Cannot change x axis scales.
		"""
		plt.figure(figsize=(self.__figure_width, self.__figure_height))

		for y in data_sets:
			plt.plot(x, y)
			if trendline:
				z = np.polyfit(x, y, 1)
				p = np.poly1d(z)
				plt.plot(x, p(x), 'r--')

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)

		if self.__y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.__y_scale, basey=self.__y_base)

		plt.show()

	def combine_data_set_names_with_x_tick_labels(self, data_set_names):
		"""
		Concatenate each data set name with each x tick label.
		"""
		new_x_tick_labels = []
		for label in self.__x_tick_labels:
			for name in data_set_names:
				new_x_tick_labels.append(label + f' ({name})')
		return new_x_tick_labels

	def multiply_x_tick_labels(self, m):
		"""
		Multiply the number of x tick labels.
		Like combine_data_set_names_with_x_tick_labels but doesn't add concatenate anything.
		"""
		new_x_tick_labels = []
		for label in self.__x_tick_labels:
			for i in range(m):
				new_x_tick_labels.append(label)
		return new_x_tick_labels

	def select_color(self, index):
		"""
		Select a color.
		"""
		return self.__COLORS[index]

	def generate_colors(self, number_of_colors, difference_threshold):
		"""
		Generate random colors.
		Filter out too similar of colors.
		"""
		n = 0
		colors = []
		while n < number_of_colors:
			new_color = self.generate_random_color()
			is_different_enough = True
			for color in colors:
				if self.color_distance(color, new_color) < difference_threshold:
					is_different_enough = False
					break
			if is_different_enough:
				colors.append(new_color)
				n += 1
		return colors

	def color_distance(self, a, b):
		"""
		Get the distance between two colors.
		Used to determine similarity.
		Note: not the best approach. 
		"""
		r = a[0] - b[0]
		g = a[1] - b[1]
		b = a[2] - b[2]
		return math.sqrt(r ** 2 + g ** 2 + b ** 2)

	def generate_random_color(self):
		"""
		Generate a random color.
		Alpha is always 1.
		"""
		return (random.random(), random.random(), random.random(), 1) #rgba values

	def generate_seed(self):
		"""
		Generate seed.
		"""
		ddate = datetime.datetime.now().date()
		year = ddate.strftime('%Y')
		month = ddate.strftime('%m')
		day = ddate.strftime('%d')
		date_str = str(month) + str(day) + str(year)
		date_int = int(date_str)
		dtime = datetime.datetime.now().time()
		hour = dtime.hour
		minute = dtime.minute
		second = dtime.second
		time_str = str(hour) + str(minute) + str(second)
		time_int = int(time_str)
		nanoseconds_since_epoch = time.time_ns()
		return  (nanoseconds_since_epoch * date_int) % time_int
		
#Example usage
if __name__ == '__main__':
	#Create some fake data
	line_plot_data_x = np.linspace(0, 5, 100)
	line_plot_data_y_a = line_plot_data_x**2
	line_plot_data_y_b = line_plot_data_x**3

	bar_plot_data_a = np.random.rand(10)
	bar_plot_data_b = np.random.rand(10)
	bar_plot_data_c = np.random.rand(10)
	bar_plot_data_d = np.random.rand(10)
	bar_plot_data_e = np.random.rand(10)

	box_plot_data_a = np.empty(10, dtype=object)
	box_plot_data_b = np.empty(10, dtype=object)
	for i in range(0, 10):
		box_plot_data_a[i] = np.random.rand(10)
		box_plot_data_b[i] = np.random.rand(20)

	histogram_data = np.random.randn(10000)

	scatter_plot_data_x = np.random.randn(10000)
	scatter_plot_data_y = np.random.randn(10000)

	heatmap_data = np.zeros((10, 10))
	for i in range(0, 10):
		for j in range(0, 10):
			heatmap_data[i][j] = i * j

	#Set up the plot generator
	pg = PlotGenerator(use_random=True)
	pg.set_figure_options(width = 10)
	pg.set_text_options(45, 'right', 0, 'center', 8)
	pg.set_text('Title', 'X', 'Y', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'], None)

	#Plot the line plot data
	#pg.line_plot(line_plot_data_x, [line_plot_data_y_a, line_plot_data_y_b], trendline=True)

	#Plot the bar plot data
	#SETS = []
	#NAMES = []
	#for i in range(1, 6):
	#	SETS.append(np.random.rand(10))
	#	NAMES.append(chr(i-1+ord('A')))
	#	pg.bar_plot(SETS, NAMES)

	#Plot the box plot data
	#Note the x tick labels. 
	#It will automatically generate this using the previously set x tick labels and given data_a_name and data_b_name values
	#pg.box_plot([box_plot_data_a, box_plot_data_b], ['A', 'B'], False)

	#Plot the histogram data
	sets = []
	for i in range(2):
		sets.append(np.random.rand(1000))

	pg.set_axis_options('linear', 10, 'log', 10)
	pg.histogram(sets, 100, ['A', 'B'])

	#Plot the scatter plot data
	#pg.scatter_plot(scatter_plot_data_x, scatter_plot_data_y, trendline=True)

	#Plot the heatmap data
	#pg.set_text('Title', 'X', 'Y', None, None)
	#pg.heatmap(heatmap_data)