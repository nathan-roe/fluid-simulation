import os
from enum import Enum
from math import floor

class FluidCell:
	def __init__(self, density=10):
		self.velocity = 0
		self.density=density
		self.x = 0
		self.y = 0


class Simulation:
	GRAVITY: float = 10.0

	class FieldType(Enum):
		U_FIELD = 0,
		V_FIELD = 1,
		S_FIELD = 2

	def __init__(self, density, over_relaxation):
		self.density=density
		self.dt = 1.0 / 120.0
		self.over_relaxation = over_relaxation
		width,height = Simulation.get_terminal_dimensions()
		self.num_x = width + 2
		self.num_y = height + 2
		self.num_cells = self.num_x * self.num_y
		self.h = height / 2
		self.u = [0.0 for _ in range(self.num_cells)]
		self.v = [0.0 for _ in range(self.num_cells)]
		self.new_u = [0.0 for _ in range(self.num_cells)]
		self.new_v = [0.0 for _ in range(self.num_cells)]
		self.p = [0.0 for _ in range(self.num_cells)]
		self.s = [0.0 for _ in range(self.num_cells)]
		self.m = [1.0 for _ in range(self.num_cells)]
		self.new_m = [0.0 for _ in range(self.num_cells)]

	def integrate(self):
		for i in range(1, self.num_x):
			for j in range(1, self.num_y - 1):
				if(self.s[i*self.num_y + j] != 0.0 and self.s[i*self.num_y + j-1] != 0.0):
					self.v[i*self.num_y + j] += Simulation.GRAVITY * self.dt



	def solve_incompressibility(self, num_iters):
		cp = self.density * self.h / self.dt
		for iter in range(num_iters):
			for i in range(1, self.num_x-1):
				for j in range(1, self.num_y - 1):
					if(self.s[i*self.num_y + j] == 0.0):
						continue
					s = self.s[i*self.num_y + j]
					sx0 = self.s[(i-1)*self.num_y + j]
					sx1 = self.s[(i+1)*self.num_y + j]
					sy0 = self.s[i*self.num_y + j-1]
					sy1 = self.s[i*self.num_y + j+1]
					s = sx0 + sx1 + sy0  + sy1
					if (s==0.0):
						continue
					div = self.u[(i+1)*self.num_y + j] - self.u[i*self.num_y + j] \
						+ self.v[i*self.num_y + j+1] - self.v[i*self.num_y + j]
					
					p = -div / s
					p *= self.over_relaxation
					self.p[i*self.num_y] += cp * p

					self.u[i *self.num_y + j] -= sx0 * p
					self.u[(i+1)*self.num_y + j] += sx1 * p
					self.v[i*self.num_y + j] -= sy0 * p
					self.v[i*self.num_y + j+1] += sy1 * p

	def extrapolate(self):
		for i in range(self.num_x):
			self.u[i*self.num_y] = self.u[i*self.num_y + 1]
			self.u[i*self.num_y + self.num_y-1] = self.u[i*self.num_y + self.num_y-2]
		
		for j in range(self.num_y):
			self.v[j] = self.v[self.num_y + j]
			self.v[(self.num_x - 1)*self.num_y + j] = self.v[(self.num_x-2)*self.num_y + j]

	def sample_field(self, x: float, y: float, field: FieldType):
		h1 = 1.0 / self.h
		h2 = 0.5 * self.h
		x = max(min(x, self.num_x * self.h), self.h)
		y = max(min(y, self.num_y * self.h), self.h)
		dx = 0.0
		dy = 0.0
		f = 0.0

		match(field):
			case Simulation.FieldType.U_FIELD:
				f = self.u
				dy = h2
			case Simulation.FieldType.V_FIELD:
				f = self.v
				dx = h2
			case Simulation.FieldType.S_FIELD:
				f = self.m
				dx = h2
		
		x0 = min(floor((x-dx)*h1), self.num_x-1)
		tx = ((x-dx) - x0*self.h) * h1
		x1 = min(x0 + 1, self.num_x-1)

		y0 = min(floor((y-dy)*h1, self.num_y-1))
		ty = ((y-dy), - y0*self.h) * h1
		y1 = min(y0 + 1, self.num_y-1)

		sx = 1.0 - tx
		sy = 1.0 - ty

		return sx*sy * f[x0*self.num_y + y0] + \
			tx*sy * f[x1*self.num_y + y0] + \
			tx*ty * f[x1*self.num_y + y1] + \
			sx*ty * f[x0*self.num_y + y1]

	def avg_u(self, i, j):
		return (self.u[(i-1)*self.num_y + j] + self.u[i*self.num_y + j] + \
			self.u[(i-1)*self.num_y + j+1] + self.u[i*self.num_y + j+1]) * 0.25
	
	def avg_v(self, i, j):
		return (self.v[(i-1)*self.num_y + j] + self.v[i*self.num_y + j] + \
			self.v[(i-1)*self.num_y + j+1] + self.v[i*self.num_y + j+1]) * 0.25

	def advect_velocity(self):
		self.new_u = self.u
		self.new_v = self.v
		h2 = 0.5 * self.h

		for i in range(1, self.num_x):
			for j in range(1, self.num_y):
				if (self.s[i*self.num_y + j] != 0.0 and self.s[i*self.num_y + j-1] != 0.0 and i < self.num_x - 1):
					x = i*self.h + h2
					y = j*self.h
					u = self.avg_u(i, j)
					v = self.v[i*self.num_y + j]
					x = x - self.dt*u
					y = y - self.dt*v
					v = self.sample_field(x,y,Simulation.FieldType.V_FIELD)
					self.new_v[i*self.num_y + j] = v

		self.u = self.new_u
		self.v = self.new_v

	def advect_smoke(self):
		self.new_m = self.m
		h2 = 0.5 * self.h

		for i in range(1, self.num_x-1):
			for j in range(1, self.num_y-1):
				if (self.s[i*self.num_y + j] != 0.0):
					u = (self.u[i*self.num_y + j] + self.u[(i+1)*self.num_y + j]) * 0.5
					v = (self.v[i*self.num_y + j] + self.v[i*self.num_y + j+1]) * 0.5
					x = i*self.h + h2 - self.dt*u
					y = j*self.h + h2 - self.dt*v
					self.new_m[i*self.num_y + j] = self.sample_field(x,y,Simulation.FieldType.S_FIELD)
		self.m = self.new_m
	
	def simulate(self, num_iters):
		self.integrate()
		self.solve_incompressibility(num_iters)
		self.extrapolate()
		self.advect_velocity()
		self.advect_smoke()

	def draw(self):
		min_p = self.p[0]
		max_p = self.p[0]

		for i in range(self.num_cells):
			min_p = min(min_p, self.p[i])
			max_p = max(max_p, self.p[i])
		l = []
		for i in range(self.num_x):
			for j in range(self.num_y):
				l.append(self.p[i*self.num_y + j])
		return l

	def get_sci_color(self):
		pass

	def get_terminal_dimensions():
		size = os.get_terminal_size()
		return size.columns, size.lines


simulation = Simulation(1000.0, 1.9)