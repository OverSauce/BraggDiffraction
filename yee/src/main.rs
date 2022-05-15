use std::io::Write;                                                                                                                                                                                                                                                                                                                              
use std::fs::File;

const EPS0: f64 = 8.854187817e-12;
const MU0: f64 = 1.2566370614e-6;
// const PI: f64 = 3.14159265358979323846;
const N: usize = 100;
const TIME: usize = 400;

fn main() {
    let mut file = File::create("output.txt").expect("create failed");

    let c: f64 = 1.0/(EPS0*MU0).sqrt();
    let mut ez: Vec<f64> = vec![0.0; N+1];
    let mut hy: Vec<f64> = vec![0.0; N];
    let x = linspace(-250.0*1e-6, 250.0*1e-6, N+1);
    let dx = x[1]-x[0];
    let dt = dx/c;
    
    let mut ez_temp: Vec<f64> = vec![0.0; N+1];
    let mut hy_temp: Vec<f64> = vec![0.0; N];

    for t in 0..TIME {
        for i in 0..N {
            hy_temp[i] = hy[i] + (ez[i+1]-ez[i])*dt/dx/MU0;
        }

        hy = hy_temp.to_vec();

        for i in 1..N {
            ez_temp[i] = ez[i] + (hy[i]-hy[i-1])*dt/dx/EPS0;
        }
        ez_temp[0] = ez[1];
        ez_temp[N] = ez[N-1]; 

        ez = ez_temp.to_vec();

        ez[N/2] = (-0.5*((t as f64 - 40.0)/10.0)*((t as f64 - 40.0)/10.0)).exp();

        for i in 0..N {
            write!(file, "{}, ", ez[i]).expect("write failed");
        }
        write!(file, "\n").expect("write failed");                                                                                                                     
    }
}

fn linspace(start: f64, end: f64, num: usize) -> Vec<f64> {
    let mut vec: Vec<f64> = Vec::new();
    let step: f64 = (end - start) / (num as f64);
    for i in 0..num {
        vec.push(start + i as f64 * step);
    }
    vec
}
