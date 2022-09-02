use fdtd::em::{EMBuilder, EM};
use inline_python::{python, Context};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let context = Context::new();

    const N: usize = 1000;
    const T: usize = 1000;

    let (x0, x1): (f64, f64) = (-0.1, 0.1);
    let eps0: f64 = 8.8541878176e-12;
    let mu0: f64 = 1.2566370614e-6;
    let vel: f64 = 1. / (eps0 * mu0).sqrt();

    let x: Vec<f64> = (0..N)
                        .into_iter()
                        .map(|i| x0 + (x1 - x0) * i as f64 / (N - 1) as f64)
                        .collect();
    let ds = x[1] - x[0];
    let dt = ds / vel;

    fn gaussian(x: f64, u: f64, sigma: f64) -> f64 {
        (-((x - u) / (sigma)).powi(2) / 2.).exp()
    }

    let mut em: EM<f64> = EMBuilder::new()
        .dimensions(T, N)
        .delta(vel, ds)
        .waveform(
            (0..N)
                .into_iter()
                .map(|s| gaussian(s as f64 * ds, N as f64 / 2. * ds, 0.01))
                .collect(),
        )
        .build();

    em.update((eps0, mu0), (dt, ds));

    #[allow(unused_variables)] // used only in python
    let ez = em.get_efield(10);
    let hy = em.get_hfield(10);

    context.run(python! {
        import numpy as np
        import matplotlib.pyplot as plt

        plt.ion()
        fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)
        for t in range(int('T / 10)):
            ax0.cla()
            ax0.plot('x, 'ez[t])
            ax0.set_ylim([-1.1, 1.1])
            ax0.set_ylabel("E field")
            ax1.cla()
            ax1.plot('x, 'hy[t])
            ax1.set_ylim([-1.5e-3, 1.5e-3])
            ax1.set_ylabel("H field")
            ax0.set_title("%.0f fs" % (1e12 * t * 'dt))
            ax1.set_xlabel("x[m]")
            plt.pause(0.001)
        plt.ioff()
    });

    drop(context); // not sure if this is needed

    Ok(())
}
