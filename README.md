# Cloud Physics: Droplet Trajectory and Microphysics

This repository implements a numerical simulation of a cloud droplet's life cycle. It models the physical behavior and trajectory of a water droplet as it moves through and below a cloud layer until it interacts with the ground surface.

## Physical Process

The simulation uses a semi-analytical approach to track a droplet's motion and mass change (growth/evaporation) throughout its lifecycle. It handles two distinct environments:
1. **Inside the Cloud**: Growth via condensation and interaction with the updraft.
2. **Below Cloud Base**: Evaporation and descent in subsaturated air.

### Key Parameters Calculated

The model provides detailed output for the following physical milestones:
- **Maximum Height**: Relative to both the starting point and the ground.
- **Microphysical Evolution**: Droplet size at maximum height, at the cloud base, and upon reaching the ground.
- **Timing metrics**: Time of ascent, time of descent, total residency time within the cloud, and total duration of the rain process.
- **Evaporation Loss**: The difference in raindrop mass/radius between the cloud base exit and ground impact.

## Visualization

The project generates several diagnostics plots to illustrate the droplet dynamics:
- `height-radius.jpg`: Relationship between vertical position and droplet size.
- `velocity-radius.jpg`: Coordination between terminal fall velocity and droplet mass.
- `height-time.jpg` & `height-velocity.jpg`: Phase-space plots of the droplet trajectory.

## Usage

Run the main simulation script:
```bash
python cloud_physics_project.py
```

## Requirements
- Python 3.x
- NumPy
- Matplotlib
