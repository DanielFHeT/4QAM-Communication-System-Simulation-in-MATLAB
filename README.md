# 4QAM Communication System Simulation in MATLAB

This repository provides a MATLAB-based simulation of a 4QAM (Quadrature Amplitude Modulation) communication system. The simulation includes both modulation and demodulation functions, signal constellation plotting, and error analysis, designed to help visualize and understand 4QAM signal processing.

## Overview

Quadrature Amplitude Modulation (QAM) is a widely used modulation scheme in digital communication systems. This project focuses on 4QAM, a modulation type that uses four distinct symbol points in the complex plane, allowing two bits to be transmitted per symbol. The simulation includes:

- **Modulation and Demodulation**: Converting bits to complex symbols for transmission and back to bits after reception.
- **Constellation Plot**: Visualization of the 4QAM symbol map, showing symbol points in the complex plane.
- **Bit Error Rate (BER) Analysis**: Calculation of the BER based on signal-to-noise ratio (SNR) to assess system performance.

## How to Run the Code in MATLAB

1. **Clone the Repository**: Clone the repository to your local machine or download it as a ZIP file.
   ```bash
   git clone https://github.com/DanielFHeT/4QAM-Communication-System-Simulation-in-MATLAB.git
Open MATLAB: Launch MATLAB and navigate to the folder containing the code files.

Run the Simulation: Run the main script main.m to start the simulation.

matlab
Copy code
main;
Provide Simulation Parameters: The script may prompt for parameters, such as the number of symbols, SNR levels, or other settings to configure the simulation according to your needs.

View the Results:

The constellation plot will display the 4QAM symbol points in the complex plane.
The BER vs. SNR plot will show how the error rate changes with varying levels of noise.
Features
Modulation and Demodulation: Implements 4QAM modulation, mapping each pair of bits to one of four complex symbols, and demodulation, converting received symbols back to bits.
Constellation Diagram: Visual representation of the 4QAM symbol constellation, aiding in understanding symbol positions.
BER Calculation: Estimates the bit error rate for different SNR values, helping to analyze system performance.
Noise Simulation: Adds AWGN (Additive White Gaussian Noise) to the transmitted signal to simulate real-world channel conditions.
Requirements
MATLAB installed on your system.
Basic knowledge of MATLAB to navigate directories and run scripts.
Understanding of digital modulation concepts is beneficial but not necessary.

License
This project is released under the MIT License, allowing free use, modification, and distribution.

Contact
For questions or suggestions, please open an issue or contact DanielFHeT.
