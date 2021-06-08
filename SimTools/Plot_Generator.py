import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#TODO: Decide what data format (list, numpy, pandas) to use with Professor Miller
class Plot_Generator:
	"""
	Class for generating plots.
	"""
	def __init__(self):
		"""
		Initialize self.
		Sets all settables to their default values.
		"""
		self.x_scale = None #matplotlib defaults to linear if not specified
		self.x_base = 10
		self.y_scale = None #matplotlib defaults to linear if not specified
		self.y_base = 10
		self.x_tick_label_rotation = 0
		self.x_tick_label_horizontal_alignment = 'center'
		self.y_tick_label_rotation = 0
		self.y_tick_label_horizontal_alignment = 'center'

		self.title = ""
		self.x_label = ""
		self.y_label = ""
		self.x_tick_labels = None
		self.y_tick_labels = None

		self.COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

	def set_text(self, title, x_label, y_label, x_tick_labels, y_tick_labels):
		"""Sets the text of plots.
		"""
		self.title = title
		self.x_label = x_label
		self.y_label = y_label
		self.x_tick_labels = x_tick_labels
		self.y_tick_labels = y_tick_labels

	def set_text_options(self, x_tick_label_rotation, x_tick_label_horizontal_alignment, y_tick_label_rotation, y_tick_label_horizontal_alignment):
		"""
		Sets the text options of plots.
		"""
		self.x_tick_label_rotation = x_tick_label_rotation
		self.x_tick_label_horizontal_alignment = x_tick_label_horizontal_alignment
		self.y_tick_label_rotation = y_tick_label_rotation
		self.y_tick_label_horizontal_alignment = y_tick_label_horizontal_alignment

	def set_axis_options(self, x_scale, x_base, y_scale, y_base):
		"""
		Set the axis options of the plots.
		Changing the bases does nothing if the scale is None or 'linear'.
		"""
		if x_scale == 'linear':
			self.x_scale = None
		else:
			self.x_scale = x_scale
		self.x_base = x_base
		if y_scale == 'linear':
			self.y_scale = None
		else:
			self.y_scale = y_scale
		self.y_base = y_base

	def bar_plot(self, data_sets, data_set_names):
		"""
		Generates a bar plot of one or many data sets.
		"""
		assert len(data_sets) == len(data_set_names)
		assert len(data_sets) <= len(self.COLORS)
		for i in range(1, len(data_sets)):
			assert len(data_sets[i]) == len(data_sets[i - 1])

		NUM_SETS = len(data_sets)

		x = np.arange(0, len(data_sets[0]))

		plt.figure()

		for i in range(0, NUM_SETS):
			self.gen_bar_plot_object(data_sets[i], self.COLORS[i] , i, NUM_SETS)

		if self.y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.y_scale, base=self.y_base)

		plt.title(self.title)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label)

		if self.x_tick_labels != None:
			plt.xticks(x, labels=self.x_tick_labels, rotation=self.x_tick_label_rotation, ha=self.x_tick_label_horizontal_alignment)

		plt.legend(data_set_names)

		plt.show()

	def gen_bar_plot_object(self, data, color, plot_num, num_plots):
		x = np.arange(0, len(data), 1)

		width = 1 / num_plots
		offset = width * plot_num - width / num_plots
		if num_plots % 2 != 0:
			offset = width * plot_num - width

		bar_plot = plt.bar(x + offset, data, color=color, width=width)
		return bar_plot



	def comparision_box_plot(self, data_a, data_b, data_a_name, data_b_name, showfliers):
		"""
		Generates a box plot comparing two data sets.
		"""
		assert len(data_a) == len(data_b)

		plt.figure()

		box_a = self.gen_box_plot_object(data_a, 'red', 0, 2, showfliers)
		box_b = self.gen_box_plot_object(data_b, 'blue', 1, 2, showfliers)

		if self.y_scale != None: #Needed because matplotlib does not like setting the base for linear plots
			plt.yscale(self.y_scale, base=self.y_base)

		plt.title(self.title)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label)
		if self.x_tick_labels != None:
			#Generate x tick labels
			new_x_tick_labels = []
			for label in self.x_tick_labels:
				new_x_tick_labels.append(label + " (" + data_a_name + ")")
				new_x_tick_labels.append(label + " (" + data_b_name + ")")

			x_ticks = np.arange(0, len(data_a) + len(data_b), 1)	
			plt.xticks(x_ticks, labels=new_x_tick_labels, rotation=self.x_tick_label_rotation, ha=self.x_tick_label_horizontal_alignment)

		plt.legend([box_a['boxes'][0], box_b['boxes'][0]], [data_a_name, data_b_name])

		plt.show()

	def gen_box_plot_object(self, data, color, plot_num, num_plots, showfliers):
		positions = np.arange(plot_num, len(data) * num_plots, num_plots)
		box_plot = plt.boxplot(data, patch_artist=True, positions=positions, showfliers=showfliers)
		for box in box_plot['boxes']:
			box.set(color=color, linewidth=1)
			box.set(facecolor='white')
		for flier in box_plot['fliers']:
			flier.set(markersize=0.5, alpha=0.5)
		return box_plot

#Example plots using Plot_Generator
if __name__ == '__main__':
	#Create some fake data
	bar_plot_data_a = []
	bar_plot_data_b = []
	for i in range(1, 11):
		bar_plot_data_a.append((i+1))
		bar_plot_data_b.append((i+1))

	box_plot_data_a = []
	box_plot_data_b = []
	for i in range(0, 10):
		box_plot_data_a.append([])
		box_plot_data_b.append([])
		for j in range(0, 10):
			box_plot_data_a[i].append((i+1) * (j+1))
			box_plot_data_b[i].append((i+1) + (j+1))

	#Set up the plot generator
	pg = Plot_Generator()
	pg.set_text_options(45, 'right', 0, 'center')
	pg.set_text('Title', 'X', 'Y', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'], None)

	#Plot the bar plot data
	pg.bar_plot([bar_plot_data_a, bar_plot_data_b], ['A', 'B'])

	#Plot the box plot data
	#Note the x tick labels. 
	#It will automatically generate this using the previously set x tick labels and given data_a_name and data_b_name values
	pg.comparision_box_plot(box_plot_data_a, box_plot_data_b, 'A', 'B', False)