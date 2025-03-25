import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import time

from constants import (SEMI_MAJOR_AXIS, ECCENTRICITY, LONGITUDE_OF_ASCENDING_NODE, 
                       ARGUMENT_OF_PERIAPSIS, INCLINATION, MEAN_ANOMALY, T, T0, J20)
from constants import PATH_FIRST_STEP, PATH_VARIATION_STEP, PATH_FIGURE, PATH_VARIATION_LOG, SUBFOLDERS
from date_helper import DateHelper
from file_helper import FileHelper
from utilities import calculate_initial_coordinate_velocity
from logger import Logger

start_time = time()

logger = Logger(PATH_VARIATION_LOG)
dh = DateHelper()
fh = FileHelper()

# Get the initial coordinates and velocities
print("Calculating initial coordinates and velocities.")
logger.log_info("Calculating initial coordinates and velocities.")
start = time()
initial_coordinate, initial_velocity = calculate_initial_coordinate_velocity(
    semi_major_axis=SEMI_MAJOR_AXIS,
    eccentricity=ECCENTRICITY,
    longitude_of_ascending_node=LONGITUDE_OF_ASCENDING_NODE,
    argument_pericenter=ARGUMENT_OF_PERIAPSIS,
    inclination=INCLINATION,
    mean_anomaly=MEAN_ANOMALY
)
end = time()
runtime = end - start
logger.log_info(f"Finished calculating intial coordinates and velicities. Elapsed time: {runtime:.6f} s.")
logger.log_info(f"Initial coordinates: {initial_coordinate}, initial velocities: {initial_velocity}")
print(f"Finished calculating intial coordinates and velicities. Elapsed time: {runtime:.6f} s.")
print(f"Initial coordinates: {initial_coordinate}, initial velocities: {initial_velocity}")

# Computing integration for initial data
logger.log_info("Calculating integration for initial data.")
print("Calculating integration for initial data. Please wait.")
t = T / 100

# Setup dataframe for saving data
df = pd.DataFrame({
    "T0": [T0],
    "t": [0],
    "x_1": [initial_coordinate[0]],
    "x_2": [initial_coordinate[1]],
    "x_3": [initial_coordinate[2]],
    "vx_1": [initial_velocity[0]],
    "vx_2": [initial_velocity[1]],
    "vx_3": [initial_velocity[2]]
})

start = time()
t_0 = dh.add_integration_time(0)
date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_0)

logger.log_info("Initial data")
logger.log_info(f"T0: {T0}")
logger.log_info(f"t: {0}")
logger.log_info(f"JD: ( {year} {month} {date} {hour} {minute} {second} )")
logger.log_info(f"Initial coordinate: {fh.get_coordinate(i=0)}")
logger.log_info(f"Initial velocity: {fh.get_velocity(i=0)}")

# Generate the initial prediction for orbit (x10, x20, x30, x'10, x'20, x'30)

t_n = dh.add_integration_time(T)
date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_n)

# Configure and run the integrator
fh.configure_file(
    coordinate=initial_coordinate,
    velocity=initial_velocity,
    date=date, 
    month=month, 
    year=year, 
    hour=hour, 
    minute=minute, 
    second=second, 
    t=t
)
fh.run_exe_file()

# Get computed coordinates and velocities
computed_coordinates, computed_velocities = fh.get_coordinates(), fh.get_velocities()

# Add the computed coordinates and velocities to df
computed_coordinates_df = pd.DataFrame(computed_coordinates, columns=["x_1", "x_2", "x_3"])
computed_velocities_df = pd.DataFrame(computed_velocities, columns=["vx_1", "vx_2", "vx_3"])
computed_df = pd.DataFrame(
    {
        "T0": [T0 for _ in range(1, 101)],
        "t": [t*i for i in range(1, 101)]
    }
)
computed_df = pd.concat([computed_df, computed_coordinates_df, computed_velocities_df], axis=1)
df = pd.concat([df, computed_df], ignore_index=True)

# Save the value to csv "first_step_data.csv"
df.to_csv(PATH_FIRST_STEP, index=False)
logger.log_info(f"Data is saved to {PATH_FIRST_STEP}\n\n")
print(f"Data is saved to {PATH_FIRST_STEP}!")
end = time()
runtime = end - start
logger.log_info(f"Finished calculating integration for intial coordinates and velicities. Elapsed time: {runtime:.6f} s.")
print(f"Finished calculating integration for intial coordinates and velicities. Elapsed time: {runtime:.6f} s.")

# Add the variation to the first component of x vector. Then generate the prediction.
logger.log_info("Calculating variations.")
print("Calculating variations. Please wait!")

start = time()
for exp in range(-12, 0):
    variation = pow(10, exp)
    for index in range(6):
        variations = np.zeros(6)
        variations[index] = variation
        variations_splitted = np.split(variations, 2)

        coordinate_with_variation = initial_coordinate + variations_splitted[0]
        velocity_with_variation = initial_velocity + variations_splitted[1]

        t_0 = dh.add_integration_time(t * 0)
        date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_0)
        
        # Setup dataframe for saving data
        df = pd.DataFrame({
            "T0": [T0],
            "t": [0],
            "x_1": [coordinate_with_variation[0]],
            "x_2": [coordinate_with_variation[1]],
            "x_3": [coordinate_with_variation[2]],
            "vx_1": [velocity_with_variation[0]],
            "vx_2": [velocity_with_variation[1]],
            "vx_3": [velocity_with_variation[2]]
        })
        logger.log_info(f"Initial data with variations {variations}")
        logger.log_info(f"T0: {T0}")
        logger.log_info(f"t: {0}")
        logger.log_info(f"JD: ( {year} {month} {date} {hour} {minute} {second} )")
        logger.log_info(f"Initial coordinate: {coordinate_with_variation}")
        logger.log_info(f"Initial velocity: {velocity_with_variation}")
        
        start_variation = time()
        
        t_n = dh.add_integration_time(T)
        date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_n)

        # Configure and run the integrator
        fh.configure_file(
            coordinate=coordinate_with_variation,
            velocity=velocity_with_variation,
            date=date, 
            month=month, 
            year=year, 
            hour=hour, 
            minute=minute, 
            second=second, 
            t=t
        )
        fh.run_exe_file()
        
        end_variation = time()
        runtime_variation = end_variation - start_variation
        logger.log_info(f"Finished calculating for variation 10^{exp} for column: {SUBFOLDERS[index]}. Elapsed time: {runtime:.6f} s.")
        print(f"Finished calculating for variation 10^{exp} for column: {SUBFOLDERS[index]}. Elapsed time: {runtime:.6f} s.")
        # Get computed coordinates and velocities
        computed_coordinates_variation, computed_velocities_variation = fh.get_coordinates(), fh.get_velocities()

        # Add the computed coordinates and velocities to df
        computed_coordinates_df = pd.DataFrame(computed_coordinates_variation, columns=["x_1", "x_2", "x_3"])
        computed_velocities_df = pd.DataFrame(computed_velocities_variation, columns=["vx_1", "vx_2", "vx_3"])
        computed_df = pd.DataFrame(
            {
                "T0": [T0 for _ in range(1, 101)],
                "t": [t*i for i in range(1, 101)]
            }
        )
        computed_df = pd.concat([computed_df, computed_coordinates_df, computed_velocities_df], axis=1)
        df = pd.concat([df, computed_df], ignore_index=True)

        # Save the value to csv "first_step_data.csv"
        path = PATH_VARIATION_STEP.format(f"{SUBFOLDERS[index]}", f"10^{exp}")
        df.to_csv(path, index=False)
        logger.log_info(f"Data is saved to {path}\n\n")
        print(f"Data is saved to {path}!")
    
    # Variation to J20
    fh.modify_garm360_in(new_value=J20 + variation)
    t_0 = dh.add_integration_time(t * 0)
    date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_0)
    
    # Setup dataframe for saving data
    df = pd.DataFrame({
        "T0": [T0],
        "t": [0],
        "x_1": [initial_coordinate[0]],
        "x_2": [initial_coordinate[1]],
        "x_3": [initial_coordinate[2]],
        "vx_1": [initial_velocity[0]],
        "vx_2": [initial_velocity[1]],
        "vx_3": [initial_velocity[2]]
    })
    logger.log_info(f"Initial data with variations on J20 {variation}")
    logger.log_info(f"T0: {T0}")
    logger.log_info(f"t: {0}")
    logger.log_info(f"JD: ( {year} {month} {date} {hour} {minute} {second} )")
    logger.log_info(f"Initial coordinate: {initial_coordinate}")
    logger.log_info(f"Initial velocity: {initial_velocity}")
    
    start_variation = time()
    
    t_n = dh.add_integration_time(T)
    date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_n)

    # Configure and run the integrator
    fh.configure_file(
        coordinate=initial_coordinate,
        velocity=initial_velocity,
        date=date, 
        month=month, 
        year=year, 
        hour=hour, 
        minute=minute, 
        second=second, 
        t=t
    )
    fh.run_exe_file()
    
    end_variation = time()
    runtime_variation = end_variation - start_variation
    logger.log_info(f"Finished calculating for variation 10^{exp} for J20. Elapsed time: {runtime:.6f} s.")
    print(f"Finished calculating for variation 10^{exp} for J20. Elapsed time: {runtime:.6f} s.")
    # Get computed coordinates and velocities
    computed_coordinates_variation, computed_velocities_variation = fh.get_coordinates(), fh.get_velocities()

    # Add the computed coordinates and velocities to df
    computed_coordinates_df = pd.DataFrame(computed_coordinates_variation, columns=["x_1", "x_2", "x_3"])
    computed_velocities_df = pd.DataFrame(computed_velocities_variation, columns=["vx_1", "vx_2", "vx_3"])
    computed_df = pd.DataFrame(
        {
            "T0": [T0 for _ in range(1, 101)],
            "t": [t*i for i in range(1, 101)]
        }
    )
    computed_df = pd.concat([computed_df, computed_coordinates_df, computed_velocities_df], axis=1)
    df = pd.concat([df, computed_df], ignore_index=True)

    # Save the value to csv "first_step_data.csv"
    path = PATH_VARIATION_STEP.format(f"{SUBFOLDERS[6]}", f"10^{exp}")
    df.to_csv(path, index=False)
    logger.log_info(f"Data is saved to {path}\n\n")
    print(f"Data is saved to {path}!")

end = time()
runtime = end - start
logger.log_info(f"Finished calculating integration for intial coordinates and velicities with variation. Elapsed time: {runtime:.6f} s.")
print(f"Finished calculating integration for intial coordinates and velicities. Elapsed time: {runtime:.6f} s.")
    
# Calculate Derivative difference (Only at the last moment of integration)
logger.log_info("Calculating Derivative difference (Only at the last moment of integration)")
print("Calculating Derivative difference (Only at the last moment of integration)")

start = time()
df = pd.read_csv(PATH_FIRST_STEP)

coordinate_velocity_columns = ["x_1", "x_2", "x_3", "vx_1", "vx_2", "vx_3"]
computed_coordinates_velocities = df[coordinate_velocity_columns].iloc[-1].to_numpy()

variations = [pow(10, exp) for exp in range(-12, 0)]
titles = [
    r"Derivation vs Variation for Coordinates $x_{1}$",
    r"Derivation vs Variation for Coordinates $x_{2}$",
    r"Derivation vs Variation for Coordinates $x_{3}$",
    r"Derivation vs Variation velocity $\dot{x_{1}}$",
    r"Derivation vs Variation velocity $\dot{x_{2}}$",
    r"Derivation vs Variation velocity $\dot{x_{3}}$"
]

for i in range(6):
    plt.clf()  # Clear the current figure
    plt.figure()  # Create a new figure

    derivatives = []
    for exp in range(-12, 0):
        var = pow(10, exp)
        path = PATH_VARIATION_STEP.format(f"{SUBFOLDERS[i]}", f"10^{exp}")
        computed_df = pd.read_csv(path)
        computed_coordinates_velocities_variation = computed_df[coordinate_velocity_columns].iloc[-1].to_numpy()

        derivative_difference = (computed_coordinates_velocities_variation - computed_coordinates_velocities) / var
        logger.log_info(f"Derivative difference for variation 10^{exp} on {SUBFOLDERS[i]}: \n{derivative_difference}")
        derivatives.append(derivative_difference)
    derivatives = np.array(derivatives)
    plt.plot(variations, derivatives[:,0], label=r"$\dfrac{\partial x_{1}}{\partial x_{10}}$")
    plt.plot(variations, derivatives[:,1], label=r"$\dfrac{\partial x_{2}}{\partial x_{20}}$")
    plt.plot(variations, derivatives[:,2], label=r"$\dfrac{\partial x_{3}}{\partial x_{30}}$")
    plt.plot(variations, derivatives[:,3], label=r"$\dfrac{\partial \dot{x_{1}}}{\partial \dot{x_{10}}}$")
    plt.plot(variations, derivatives[:,4], label=r"$\dfrac{\partial \dot{x_{2}}}{\partial \dot{x_{20}}}$")
    plt.plot(variations, derivatives[:,5], label=r"$\dfrac{\partial \dot{x_{3}}}{\partial \dot{x_{20}}}$")
    plt.xscale("log")  # Log scale for variation
    plt.ylabel("Derivation")
    plt.xlabel("Variation")
    plt.title(titles[i])
    plt.legend()
    plt.grid(True)
    # Set x-axis ticks for all data points
    plt.xticks(variations, rotation=45)
    plt.savefig(PATH_FIGURE.format(f"{SUBFOLDERS[i]}"))
    logger.log_info(f"Plot saved to {PATH_FIGURE.format(SUBFOLDERS[i])}")
    
# Derivatives for J20
plt.clf()  # Clear the current figure
plt.figure()  # Create a new figure
derivatives = []
for exp in range(-12, 0):
    var = pow(10, exp)
    path = PATH_VARIATION_STEP.format(f"{SUBFOLDERS[6]}", f"10^{exp}")
    computed_df = pd.read_csv(path)
    computed_coordinates_velocities_variation = computed_df[coordinate_velocity_columns].iloc[-1].to_numpy()

    derivative_difference = (computed_coordinates_velocities_variation - computed_coordinates_velocities) / var
    logger.log_info(f"Derivative difference for variation 10^{exp} on {SUBFOLDERS[6]}: \n{derivative_difference}")
    derivatives.append(derivative_difference)
derivatives = np.array(derivatives)
plt.plot(variations, derivatives[:,0], label=r"$\dfrac{\partial x_{1}}{\partial J_{20}}$")
plt.plot(variations, derivatives[:,1], label=r"$\dfrac{\partial x_{2}}{\partial J_{20}}$")
plt.plot(variations, derivatives[:,2], label=r"$\dfrac{\partial x_{3}}{\partial J_{20}}$")
plt.plot(variations, derivatives[:,3], label=r"$\dfrac{\partial \dot{x_{1}}}{\partial J_{20}}$")
plt.plot(variations, derivatives[:,4], label=r"$\dfrac{\partial \dot{x_{2}}}{\partial J_{20}}$")
plt.plot(variations, derivatives[:,5], label=r"$\dfrac{\partial \dot{x_{2}}}{\partial J_{20}}$")
plt.xscale("log")  # Log scale for variation
plt.ylabel("Derivation")
plt.xlabel("Variation")
plt.title(r"Derivation vs Variation $J_{20}$")
plt.legend()
plt.grid(True)
# Set x-axis ticks for all data points
plt.xticks(variations, rotation=45)
plt.savefig(PATH_FIGURE.format(f"{SUBFOLDERS[6]}"))

end = time()
runtime = end - start
logger.log_info(f"Finished calculating derivative difference. Elapsed time: {runtime:.6f} s.")
print(f"Finished calculating derivative difference. Elapsed time: {runtime:.6f} s.")

end_time = time()
elapsed_time = end_time - start_time
logger.log_info(f"Finished running variation_selection.py. Elapsed time: {elapsed_time:.6f} s.")
print(f"Finished running variation_selection.py. Elapsed time: {elapsed_time:.6f} s.")