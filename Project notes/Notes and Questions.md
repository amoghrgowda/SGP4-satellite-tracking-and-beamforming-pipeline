# questions
 - Azimuth and elevation angles??
# Intro 
- An antenna element is the smallest individual antenna that can receive and transmit electromagnetic signals.
- These elements are arranged in a grid. Each of these elements pick up signals (with small delay between each other due to slight increase in travel distance of the wave in an angle). This would mean that one element's peak might get cancelled out by another element's trough. Therefore, the receiver device applies a small delay to each element so that every element's peak and trough match exactly with other element's peak and trough. This will amplify the signal and increase the signal-to-noise ratio (SNR) dramatically.
## Advantages of an array of antenna elements (phased array)
- We can electronically steer the beam (without having to physically steer the antenna) by simply increasing or decreasing the delay.
- We can also send out multiple beams from the same array containing hundreds of elements.
- Also, we CANNOT have a huge/single antenna/element that is big because it would dissipate the signal. We need the antenna to be about the size of half of the wavelength of the wave. This is simply how physics work. Different sections of this massive sheet would end up with opposite voltages and the signal emitted by these opposing voltages cancel each other out.
>[!note] Why use array antenna instead of monopole or a direction-oriented antenna like the yagi uda?
>Electronic steering instead of a physical motor-steering

>[!note]
>Waves emitted peak in the direction perpendicular to the flow of current. Zero emitted parallel to it.

# phase 1
## TLE
- TLE (Two line element) is a fixed-width 69 character per line ASCII packet.
- It has the following components: Metadata (like the satellite ID, launch info, etc), Time Anchor (in epoch), Orbital shape and orientation, and decay modelling (caused by atmospheric drag).
## SGP4
- We could use Runge-Kutta/cowell's method to identify the satellite's position at any given time (t + delta t) by calculating the gravitational force, atmospheric drag, solar radiation pressure, third-body gravity of sun/moon, etc.
- But, this is computationally very expensive and slow, although it is a very precise calculation.
- To track multiple satellite positions, we use a lightweight O(1) constant time position calculation using SGP4. Can track thousands of satellite simultaneously.
- SGP4 is closely tied with TLE (two line element) format.
- SGP4 basically translates the TLE input data into an output structure (StateVector) across time.
- The mathematical models expect data to be in canonical units like angles in radians, and mean motion in radians/minute. Therefore, we need to make sure all data is consistent.
## Code explanation
### SGP4Types.hpp:
- Acts as a central domain model and data contract for the entire satellite tracking engine.
- Contains no execution logic or algorithms. Only defines the data structures to represent the 3D spatial points.
- Contains the constants as defined by WGS72 (for earth's radius, gravitation, etc).