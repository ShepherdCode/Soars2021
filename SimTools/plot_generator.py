import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import numpy as np
import numbers
import random

class PlotGenerator:
	"""
	Class for generating plots using matplotlib.
	"""
	def __init__(self, starting_color_index=None):
		"""
		Initialize self.
		Sets all settables to their default values.
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

		self.__COLORS = []
		for name in mcd.XKCD_COLORS:
			self.__COLORS.append(name)
		if starting_color_index == None:
			self.__STARTING_COLOR_INDEX = random.randint(0, len(self.__COLORS))
		else:
			self.__STARTING_COLOR_INDEX = starting_color_index
		#self.__COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

	def set_text(self, title, x_label, y_label, x_tick_labels, y_tick_labels):
		"""
		Sets the titles and labels of to-be-generated plots.
		"""
		if x_tick_labels != None:
			assert isinstance(x_tick_labels, list)
		if y_tick_labels != None:
			assert isinstance(y_tick_labels, list)
		self.__title = title
		self.__x_label = x_label
		self.__y_label = y_label
		self.__x_tick_labels = x_tick_labels
		self.__y_tick_labels = y_tick_labels

	def set_text_options(self, x_tick_label_rotation, x_tick_label_horizontal_alignment, y_tick_label_rotation, y_tick_label_horizontal_alignment):
		"""
		Sets the title and label text options of to-be-generated plots.
		"""
		self.__x_tick_label_rotation = x_tick_label_rotation
		self.__x_tick_label_horizontal_alignment = x_tick_label_horizontal_alignment
		self.__y_tick_label_rotation = y_tick_label_rotation
		self.__y_tick_label_horizontal_alignment = y_tick_label_horizontal_alignment

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
		assert isinstance(data_sets, (list, np.ndarray))
		for data_set in data_sets:
			assert isinstance(data_set, np.ndarray)
		NUM_SETS = len(data_sets)
		assert isinstance(data_set_names, list)
		assert NUM_SETS == len(data_set_names)
		assert NUM_SETS <= len(self.__COLORS)
		for i in range(1, NUM_SETS):
			assert len(data_sets[i]) == len(data_sets[i - 1])
		for ds in data_sets:
			assert len(ds) == len(self.__x_tick_labels)

		plt.figure()

		for i in range(0, NUM_SETS):
			self.__gen_bar_plot_object(data_sets[i], self.select_color(i) , i, NUM_SETS)

		if self.__y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.__y_scale, basey=self.__y_base)

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)
		if self.__x_tick_labels != None:
			new_x_tick_labels = self.combine_data_set_names_with_x_tick_labels(data_set_names)
			x_ticks = np.arange(0, NUM_SETS * len(data_sets[0]), 1)
			plt.xticks(x_ticks, labels=new_x_tick_labels, rotation=self.__x_tick_label_rotation, ha=self.__x_tick_label_horizontal_alignment)

		plt.legend(data_set_names, loc='upper left')

		plt.show()

	def __gen_bar_plot_object(self, data, color, plot_num, num_plots):
		"""
		Used in bar_plot function.
		"""
		#Note: assertions all handled earlier in bar_plot
		positions = np.arange(plot_num, len(data) * num_plots, num_plots)
		bar_plot = plt.bar(positions, data, color=color)
		return bar_plot

	def box_plot(self, data_sets, data_set_names, showfliers):
		"""
		Generates a box plot of one or many data sets.
		Cannot change x axis scale.
		"""
		assert isinstance(data_sets, (list, np.ndarray))
		for data_set in data_sets:
			assert isinstance(data_set, np.ndarray)
		NUM_SETS = len(data_sets)
		assert isinstance(data_set_names, list)
		assert NUM_SETS == len(data_set_names)
		assert NUM_SETS <= len(self.__COLORS)
		for i in range(1, NUM_SETS):
			assert len(data_sets[i]) == len(data_sets[i - 1])
		for ds in data_sets:
			assert len(ds) == len(self.__x_tick_labels)

		plt.figure()

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
			plt.xticks(x_ticks, labels=new_x_tick_labels, rotation=self.__x_tick_label_rotation, ha=self.__x_tick_label_horizontal_alignment)

		for_legend = []
		for box in boxes:
			for_legend.append(box['boxes'][0])
		plt.legend(for_legend, data_set_names, loc='upper left')

		plt.show()

	def __gen_box_plot_object(self, data, color, plot_num, num_plots, showfliers):
		"""
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

	def histogram(self, data, bins):
		"""
		Generates a histogram.
		Cannot change x axis scale.
		"""
		assert isinstance(data, np.ndarray)
		plt.figure()

		plt.hist(data, bins)

		if self.__y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.__y_scale, basey=self.__y_base)

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)

		plt.show()

	def scatter_plot(self, data_x, data_y, trendline=False):
		"""
		Generates a scatter plot.
		Cannot change x or y axis scale.
		"""
		assert isinstance(data_x, np.ndarray)
		assert isinstance(data_y, np.ndarray)

		plt.figure()

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
		assert isinstance(data_sets, (list, np.ndarray))
		for data_set in data_sets:
			assert isinstance(data_set, np.ndarray)
		NUM_SETS = len(data_sets)
		if NUM_SETS > 1:
			for i in range(1, NUM_SETS):
				assert len(data_sets[i]) == len(data_sets[i - 1])

		plt.figure()

		plt.imshow(data_sets, cmap='hot', interpolation='nearest')

		plt.title(self.__title)
		plt.xlabel(self.__x_label)
		plt.ylabel(self.__y_label)

		if self.__x_tick_labels != None:
			assert len(self.__x_tick_labels) == len(data_sets[0])
			plt.xticks(np.arange(len(data_sets[0])), self.__x_tick_labels, rotation=self.__x_tick_label_rotation, ha=self.__x_tick_label_horizontal_alignment)
		if self.__y_tick_labels != None:
			assert len(self.__y_tick_labels) == NUM_SETS
			plt.yticks(np.arange(len(self.__y_tick_labels)), self.__y_tick_labels, rotation=self.__y_tick_label_rotation, ha=self.__y_tick_label_horizontal_alignment)

		plt.show()

	def line_plot(self, x, data_sets, trendline=False):
		"""
		Generates a line plot.
		Cannot change x axis scales.
		"""
		plt.figure()

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
		Add the names of data sets to each x tick label.
		"""
		new_x_tick_labels = []
		for label in self.__x_tick_labels:
			for name in data_set_names:
				new_x_tick_labels.append(label + f'({name})')
		return new_x_tick_labels

	def select_color(self, index):
		"""
		Select a xkcd color. 
		"""
		selection_index = index + self.__STARTING_COLOR_INDEX
		if selection_index >= len(self.__COLORS):
			selection_index -= len(self.__COLORS)
		return self.__COLORS[selection_index]

#Example plots using PlotGenerator
if __name__ == '__main__':
	#Create some fake data
	line_plot_data_x = np.linspace(0, 5, 100)
	line_plot_data_y_a = line_plot_data_x**2
	line_plot_data_y_b = line_plot_data_x**3

	bar_plot_data_a = np.zeros(10)
	bar_plot_data_b = np.zeros(10)
	for i in range(0, 10):
		bar_plot_data_a[i] = i
		bar_plot_data_b[i] = i**2

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
	pg = PlotGenerator()
	pg.set_text_options(45, 'right', 0, 'center')
	pg.set_text('Title', 'X', 'Y', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'], None)

	#Plot the line plot data
	pg.line_plot(line_plot_data_x, [line_plot_data_y_a, line_plot_data_y_b], trendline=True)

	#Plot the bar plot data
	pg.bar_plot([bar_plot_data_a, bar_plot_data_b], ['A', 'B'])

	#Plot the box plot data
	#Note the x tick labels. 
	#It will automatically generate this using the previously set x tick labels and given data_a_name and data_b_name values
	pg.box_plot([box_plot_data_a, box_plot_data_b], ['A', 'B'], False)

	#Plot the histogram data
	pg.set_axis_options('log', 10, 'log', 10)
	pg.histogram(histogram_data, 20)

	#Plot the scatter plot data
	pg.scatter_plot(scatter_plot_data_x, scatter_plot_data_y, trendline=True)

	#Plot the heatmap data
	pg.set_text('Title', 'X', 'Y', None, None)
	pg.heatmap(heatmap_data)