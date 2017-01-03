import math
import pyglet
import random
import sys
import xml.etree.ElementTree as xmltree



class Point:
	components = ('x', 'y')
	dimension = 2
	eps = 1e-6

	def __init__(self, *args, **kwargs):
		data = [0.0] * 2
		for i in range(min(len(args), self.dimension)):
			data[i] = args[i]
		for c, v in kwargs.items():
			data[components.index(c)] = v
		self._data = tuple(data)

	@classmethod
	def close(Self, lhs, rhs):
		return abs(lhs - rhs) < Self.eps

	def normalisation(self, onZero = None):
		norm = self.size()
		if self.close(norm, 0.0):
			if onZero is not None:
				normalised = onZero
			else:
				raise ValueError("Zero-length vector cannot be normalised.")
		else:
			normalised = self / norm
		return normalised, norm

	def size(self):
		return math.sqrt(self.size2())

	def size2(self):
		return sum([v ** 2 for v in self._data])

	def __add__(self, rhs):
		return Point(*map(lambda l, r: l + r, self._data, rhs._data))

	def __mul__(self, rhs):
		return Point(*map(lambda v: v * rhs, self._data))

	def __sub__(self, rhs):
		return Point(*map(lambda l, r: l - r, self._data, rhs._data))

	def __truediv__(self, rhs):
		return Point(*map(lambda v: v / rhs, self._data))




# class Point:
	# def __init__(self, x = 0.0, y = 0.0):
		# self.x = x
		# self.y = y

	# def __getitem__(self, index):
		# return getattr(self, self._attrName(index))

	# def __setitem__(self, index, value):
		# return setattr(self, self._attrName(index), value)

	# def __add__(lhs, rhs):
		# return Point(
			# lhs.x + rhs.x,
			# lhs.y + rhs.y
		# )

	# def __sub__(lhs, rhs):
		# return Point(
			# lhs.x - rhs.x,
			# lhs.y - rhs.y
		# )

	# @staticmethod
	# def _attrName(index);
		# return {
		 # 0: 'x', 'x': 'x',
		 # 1: 'y', 'y': 'y'
	 # }[index]



class PlacerConfig:
	class Village:
		def __init__(self, node):
			self.count = int(node.get('count', 1))
			self.force = float(node.find('force').text)

	def __init__(self, filePath):
		root = xmltree.parse(filePath).getroot()
		area = root.find('area')
		self.areaRadius = float(area.find('radius').text)
		self.villages = [Village(village) for village in root.findall('village')]



class Placer:
	class Village:
		def __init__(self, config):
			self.force = config.force
			self.position = None

	def __init__(self, config):
		self.areaRadius = config.areaRadius
		self.villages = []
		for v in config.villages:
			self.villages.extend(createVillages(v))
		for v in self.villages:
			self.initialiseVillagePosition(v)

	def run(self):
		window = pyglet.window.Window(
			caption = "WorldKit:Placer",
			width = 512,
			height = 512,
			resizable = True
		)
		pyglet.clock.schedule(self.step, 1/60.0)
		pyglet.app.run()

	def step(self):
		forces = {v : self.computeTotalForce(v) for v in self.villages}

	def createVillages(self, config):
		return [Village(config) for i in range(config.count)]

	def initialiseVillagePosition(self, village):
		if village.position is not None:
			return
		while True:
			position = Point(
				random.uniform(-self.areaRadius, self.areaRadius),
				random.uniform(-self.areaRadius, self.areaRadius)
			)
			if self.areaContains(position):
				break
		village.position = position

	def areaContains(self, point):
		return point.size2() <= self.areaRadius ** 2

	def computeTotalForce(self, subject):
		force = Point()
		for partner in self.villages:
			if partner is subject:
				continue
			force += self.computeForce(source = partner, target = subject)
		return force

	def computeForce(self, source, target):
		difference = target.position - source.position
		difference, size = difference.normalisation(keepZero = True)
		return difference * (self.areaRadius * 2 - size) * source.force



def printHelp(out = sys.stdout):
	pass



def main(argv):
	try:
		arg = argv[1]
	except:
		arg = 'placement.xml'
	if arg in {'-h', '--help'}:
		printHelp()
		return
	Placer(PlacerConfig(argv[1])).run()
