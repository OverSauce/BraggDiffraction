pub mod em {

	pub struct EM<T> {
		hy: Vec<Vec<T>>,
		ez: Vec<Vec<T>>,
		times: usize,
		spaces: usize,
		eps: f64,
		mu: f64,
	}

	pub struct Simulate {
		times: usize,
		spaces: usize,
		eps: f64,
		mu: f64,
		dt: f64,
		ds: f64,
		waveform: Vec<f64>,
	}

	#[allow(dead_code)]
	impl EM<f64> {
		pub fn new(tlike: usize, slike: usize) -> EM<f64> {
			EM {
				hy: vec![vec![0.0; slike]; tlike],
				ez: vec![vec![0.0; slike]; tlike],
				times: tlike,
				spaces: slike,
				eps: 8.8541878176e-12,
				mu: 1.2566370614e-6,
			}
		}
		fn copy(&self) -> EM<f64> {
			EM {
				hy: self.hy.clone(),
				ez: self.ez.clone(),
				times: self.times,
				spaces: self.spaces,
				eps: self.eps,
				mu: self.mu,
			}
		}
		pub fn update(&mut self, (dt, ds): (f64, f64)) {
			for t in 1..self.hy.len() {
				for s in 0..self.hy[t].len() - 1 {
					self.hy[t][s] = self.hy[t - 1][s]
						+ (self.ez[t - 1][s + 1] - self.ez[t - 1][s]) * dt / ds / self.mu;
				}
				self.hy[t][self.spaces - 1] = self.hy[t - 1][self.spaces - 2];
				for s in 1..self.hy[t].len() {
					self.ez[t][s] =
						self.ez[t - 1][s] + (self.hy[t][s] - self.hy[t][s - 1]) * dt / ds / self.eps;
				}
				self.ez[t][0] = self.ez[t - 1][1];
			}
		}
		pub fn initial_cond(&mut self, waveform: &Vec<f64>) {
				self.ez[0] = waveform.clone();
		}
		pub fn get_efield(&self, everyn: usize) -> Vec<Vec<f64>> {
				(0..self.times / everyn)
						.map(|t| self.ez[t * everyn].clone())
						.collect()
		}
		pub fn get_hfield(&self, everyn: usize) -> Vec<Vec<f64>> {
				(0..self.times / everyn)
						.map(|t| self.hy[t * everyn].clone())
						.collect()
		}
	}

	impl Simulate {
		pub fn new() -> Simulate {
			Simulate {
				times: 0,
				spaces: 0,
				eps: 8.8541878176e-12,
				mu: 1.2566370614e-6,
				dt: 0.0,
				ds: 0.0,
				waveform: vec![],
			}
		}
		pub fn dimensions(mut self, times: usize, spaces: usize) -> Simulate {
			self.times = times;
			self.spaces = spaces;
			self
		}
		pub fn consts(mut self, eps: f64, mu: f64) -> Simulate {
			self.eps = eps;
			self.mu = mu;
			self
		}
		pub fn dt(mut self, dt: f64) -> Simulate {
			self.dt = dt;
			self
		}
		pub fn ds(mut self, ds: f64) -> Simulate {
			self.ds = ds;
			self
		}
		pub fn delta(mut self, vel: f64, ds: f64) -> Simulate {
			self.dt = ds / vel;
			self.ds = ds;
			self
		}
		pub fn waveform(mut self, waveform: Vec<f64>) -> Simulate {
			self.waveform = waveform;
			self
		}
		pub fn build(self) -> EM<f64> {
			let mut em = EM::new(self.times, self.spaces);
			em.initial_cond(&self.waveform);
			em.update((self.dt, self.ds));
			em
		}
	}
}