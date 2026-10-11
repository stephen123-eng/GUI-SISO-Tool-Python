# Python SISO Tool (Control System Designer)

A Python-based graphical control system design and analysis tool built using `tkinter`, `matplotlib`, and the `python-control` library. 

This project is a Python port of the graphical SISO tool originally developed by [Erivelton Gualter](https://github.com/EriveltonGualter/octave-control) for GNU Octave.

---

## Features

- **Interactive Root Locus:** Add, delete, or move poles and zeros dynamically.
- **Frequency Domain Analysis:** Full Bode diagrams (Magnitude & Phase with stability margins) and Nyquist plots.
- **Multiple System Architectures:** Support for 6 distinct control loop configurations complete with visual architecture preview diagrams.
- **Design Management:** 
  - Save multiple controller designs.
  - Toggle visibility to compare response curves.
  - **"Make Active"** feature to load any saved design back into the active workspace for further tuning.

---

## Prerequisites & Dependencies

Make sure you have Python installed along with the required libraries:

```bash
pip install numpy matplotlib control
