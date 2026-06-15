# 2D Pressure Diffusivity Simulation

MATLAB coursework code for pressure diffusivity in a homogeneous 2D porous medium.

The program implements:

- Explicit FTCS finite-difference scheme
- Implicit backward-Euler finite-difference scheme
- User inputs for domain size, grid resolution, time step, and final time
- Stability check for the explicit method
- 2D pressure contour plots at multiple time steps
- A-A centre-line pressure profiles
- Explicit vs implicit runtime and profile comparison
- Grid convergence study using the implicit solver

## Coursework Boundary Conditions

The domain is a rectangle of length `L` and width `W`.

- Top-left half: known pressure inlet, `Pin = 1e5 Pa`
- Bottom-right half: known pressure outlet, `Pout = 0 Pa`
- Left, right, top-right half, and bottom-left half: no-flow boundaries

The coursework figure legend says the green boundary is a known flow-rate boundary, but the caption states `Pin = 1e5 Pa`. This implementation follows the caption and treats the green boundary as a known-pressure inlet.

## How To Run

Open MATLAB in this folder and run:

```matlab
numerical_simulation_pressure_diffusivity
```

Press Enter at each prompt to use the baseline defaults:

- `L = 100 m`
- `W = 50 m`
- `20 x 10` grid cells
- Stable time step chosen from the explicit FTCS stability limit
- `t_final = 10 s`

Generated plots are saved to the `results/` folder.

## Main Equations

The governing equation is:

```text
dP/dt = alpha * (d2P/dx2 + d2P/dy2)
alpha = k / (phi * mu * c)
```

The explicit scheme uses:

```text
P(i,j)^(n+1) = P(i,j)^n
             + rx * (P(i+1,j)^n - 2P(i,j)^n + P(i-1,j)^n)
             + ry * (P(i,j+1)^n - 2P(i,j)^n + P(i,j-1)^n)
```

where:

```text
rx = alpha * dt / dx^2
ry = alpha * dt / dy^2
```

The explicit stability condition is:

```text
rx + ry <= 0.5
```

The implicit method solves the backward-Euler linear system:

```text
(1 + 2rx + 2ry) P(i,j)^(n+1)
- rx * P(i+1,j)^(n+1)
- rx * P(i-1,j)^(n+1)
- ry * P(i,j+1)^(n+1)
- ry * P(i,j-1)^(n+1)
= P(i,j)^n
```

## Before Coursework Submission

Add both students' full names and student numbers at the top of the MATLAB file and in the report.

## Suggested GitHub Commands

```bash
git init
git add .
git commit -m "Add pressure diffusivity simulation"
git branch -M main
git remote add origin <your-github-repo-url>
git push -u origin main
```
