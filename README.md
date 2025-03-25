# Orbit-Improvements

Labortatory assignment to improve zonal harmonic distortion coefficient.

## Files description

1. Python Files:

    * `constants.py`: Defines various constants used throughout the project.
    * `date_helper.py`: Contains helper functions for date manipulations.
    * `file_helper.py`: Provides functions to handle file operations.
    * `logger.py`: Implements logging functionality to track events and errors.
    * `main.py`: The main script that orchestrates the execution of the project.
    * `utilities.py`: Offers additional utility functions to support the main operations.
    * `variation_selection.py`: Handles the selection of variations in the orbit improvement process.

2. Text Files:

    * `requirements.txt`: Lists the Python dependencies required to run the project.

3. Other files:  

    * `EPH.OUT`: Output file from integrator.
    * `iszm_puc.in`: Input file for integrator.
    * `iszm_puc.exe`: Executable file to run the integrator.
    * `MODTOOLS/garm360.in`: Input file for zonal harmonic coefficient.

4. Folders:

    * `data`: All saved data points results.
        * `data/orbit_improvements`: All saved data points results from the `main.py` code.
        * `data/variation_selection`: All saved data points results from the `variation_selection.py` code.
    * `pics`: Graphs to help determine the correct variation.
    * `log`: All the `.log` files, to record all actions from the program.
        * `orbit_improvements.log`: Log file where all recorded actions from running the `main.py` file.
        * `variation_selection.log`: Log file where all recorded actions from running the `variation_selection.py` file.

## How to run?

1. Make sure all necessary libraries are installed. Run this command on terminal:

```sh
pip install -r requirements.txt
```

2. You can change your satelite orbital parameter inside the file `constants.py`

```python
# Orbital elements for Satellite
SEMI_MAJOR_AXIS = "SATELITE'S SEMI MAJOR AXIS"
ECCENTRICITY = "SATELITE'S ECCENTRICITY"
INCLINATION = "SATELITE'S INCLINATION"
LONGITUDE_OF_ASCENDING_NODE = "SATELITE'S LONGITUDE OF ASCENDING NODE"
ARGUMENT_OF_PERIAPSIS = "SATELITE'S ARGUMENT OF PERIAPSIS"
MEAN_ANOMALY = "SATELITE'S MEAN ANOMALY"
T = "SATELITE'S ORBITAL PERIOD"
N = "NUMBER OF INTEGRATION"
```

3. Run the file `variation_selection.py` first. Run this command:

```sh
python variation_selection.py
```

4. Check on folder `pics`, for the image results. Pick the right value of variation, then change the value in the file `main.py`.

```python
class OrbitalImprovement:
    def __init__(self):
        # Constants paremeter
        self.variation_coordinate = "VARIATION FOR COORDINATE"
        self.variation_velocity = "VARIATION FOR VELOCITY"
        self.variation_j20 = "VARIATION FOR J20"
        # Change also in the diagonal elements, the first 3 should be the same as variation_coordinate, and the other 3 would be the variation_velocity. The last row should all be 0e0. Here is the example:
        self.variations = np.array([
            [10e-4, 0e0, 0e0, 0e0, 0e0, 0e0],
            [0e0, 10e-4, 0e0, 0e0, 0e0, 0e0],
            [0e0, 0e0, 10e-4, 0e0, 0e0, 0e0],
            [0e0, 0e0, 0e0, 10e-6, 0e0, 0e0],
            [0e0, 0e0, 0e0, 0e0, 10e-6, 0e0],
            [0e0, 0e0, 0e0, 0e0, 0e0, 10e-6]
        ])
```

5. Run the `main.py` file, by running this command on the terminal:

```sh
python main.py
```

6. Your results is located at `data/orbit_improvements/result.csv`, in the terminal you will also see on which iteration the process is stopped, because the mean square of weight error is not decreasing anymore.