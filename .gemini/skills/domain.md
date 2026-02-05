# Cloud Physics Domain Rules

This file defines domain constraints and validation heuristics for cloud droplet
trajectory, growth, and evaporation models used in this project.

## Physical Assumptions
- Treat the model as a simplified bulk microphysics trajectory model.
- Use the same parameter constants as the simulation unless explicitly changed.
- Maintain non-negative droplet radius and physically plausible velocities.

## Invariants and Sanity Checks
- Radius must remain non-negative at all times.
- Height should increase during ascent, decrease during descent and fall.
- Time must be strictly increasing.
- Velocity should be non-negative.

## Manufactured-Solution Guidance (MMS)
- If adding or refactoring equations, prefer validating with a simple analytic
  case (e.g., constant velocity or linear radius change) to verify integration.

## Noether / Conservation Notes
- This simplified model does not enforce strict conservation laws, but changes
  should not introduce unbounded energy or mass growth without physical basis.
