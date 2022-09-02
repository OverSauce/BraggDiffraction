pub mod em {

	pub struct EM<T> {
		hy: Vec<Vec<T>>,
		ez: Vec<Vec<T>>,
		times: usize,
		spaces: usize,
	}

	pub struct EMBuilder {
		times: usize,
		spaces: usize,
		eps0: f64,
		mu0: f64,
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
			}
		}
		fn copy(&self) -> EM<f64> {
			EM {
				hy: self.hy.clone(),
				ez: self.ez.clone(),
				times: self.times,
				spaces: self.spaces,
			}
		}
		pub fn update(&mut self, (eps0, mu0): (f64, f64), (dt, ds): (f64, f64)) {
			for t in 1..self.hy.len() {
				for s in 0..self.hy[t].len() - 1 {
					self.hy[t][s] = self.hy[t - 1][s]
						+ (self.ez[t - 1][s + 1] - self.ez[t - 1][s]) * dt / ds / mu0;
				}
				self.hy[t][self.spaces - 1] = self.hy[t - 1][self.spaces - 2];
				for s in 1..self.hy[t].len() {
					self.ez[t][s] =
						self.ez[t - 1][s] + (self.hy[t][s] - self.hy[t][s - 1]) * dt / ds / eps0;
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

	impl EMBuilder {
		pub fn new() -> EMBuilder {
			EMBuilder {
				times: 0,
				spaces: 0,
				eps0: 8.8541878176e-12,
				mu0: 1.2566370614e-6,
				dt: 0.0,
				ds: 0.0,
				waveform: vec![],
			}
		}
		pub fn dimensions(mut self, times: usize, spaces: usize) -> EMBuilder {
			self.times = times;
			self.spaces = spaces;
			self
		}
		pub fn eps0(mut self, eps0: f64) -> EMBuilder {
			self.eps0 = eps0;
			self
		}
		pub fn mu0(mut self, mu0: f64) -> EMBuilder {
			self.mu0 = mu0;
			self
		}
		pub fn dt(mut self, dt: f64) -> EMBuilder {
			self.dt = dt;
			self
		}
		pub fn ds(mut self, ds: f64) -> EMBuilder {
			self.ds = ds;
			self
		}
		pub fn delta(mut self, vel: f64, ds: f64) -> EMBuilder {
			self.dt = ds / vel;
			self.ds = ds;
			self
		}
		pub fn waveform(mut self, waveform: Vec<f64>) -> EMBuilder {
			self.waveform = waveform;
			self
		}
		pub fn build(self) -> EM<f64> {
			let mut em = EM::new(self.times, self.spaces);
			em.initial_cond(&self.waveform);
			em.update((self.eps0, self.mu0), (self.dt, self.ds));
			em
		}
	}
}