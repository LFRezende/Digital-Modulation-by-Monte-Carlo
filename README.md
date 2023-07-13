# Digital Modulation by Monte Carlo's Algorithm - Maximum Spectral Density within Probability Threshold.

This project is the resolution of the Final Exam of the discipline of EET-50: Communication Systems, of the
Aeronautics Institute of Technology (ITA), in Brazil.

The original Final Exam was graded with a perfect score of 10/10 by the professor of EET-50.


## Introduction

In digital communications, it's important to analyze the limitations of your channel whilst guaranteeing a bare minimum of reliability of the information being transmitted.

In that sense, this project:

- Discusses the experimental modulation and their theoretical predictions for SER (Signal Error Rate) at the same plot of Eb/N0 (Energy of Bit vs Noise Energy), doing this for 4 different modulation techniques: 2-PAM, 4-PAM, 4-QAM, and 2-FSK, to prove the source code for modulation is responsive - this is the code of the Q1 file.
- Compares the signal's BER (Bit Error Rate) for Natural Mapping and Gray Mapping for the constellations for 4-QAM (Q2 file).

With the introduction finished, the core project begins:

- At the 3rd File, Q3, the algorithmic modulation must fit criteria of 1% Signal Error Rate (at maximum), whilst always choosing the most robust modulation technique available for the current ratio of Bit Energy and Noise Energy (Eb/N0). The technique evaluated are BPSK, QPSK, 8PSK, 16-QAM, 64-QAM. This project envisions to maximize the spectral efficiency of the data (optimizing the rate of data transmission).
- At the 4th File, Q4, basically the same, but now we study the Bit Error Rate instead of SER. 

I published this to hopefully help someone with the same curriculum required.

If it is helpful to you, let me know about it! https://twitter.com/lf_rezende_py

# Q1 - SER comparison in Experimental Digital Modulation and its Theoretical Analog

# Q2 - Bit Error Rate comparison for QAM modulation - Gray vs Natural Mapping

# Q3 - Maximizing Spectral Density with threshold of Maximum Probability of Error

# Q4 - Maximizing Spectral Density with threshold of Maximum Probability of Error
