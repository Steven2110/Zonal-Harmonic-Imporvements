import pandas as pd
import numpy as np
from time import time

from logger import Logger
from constants import (SEMI_MAJOR_AXIS, ECCENTRICITY, LONGITUDE_OF_ASCENDING_NODE, 
                       ARGUMENT_OF_PERIAPSIS, INCLINATION, MEAN_ANOMALY, T, T0, N)
from constants import (ARCSECOND, J20)
from constants import PATH_ORBIT_LOG, PATH_OBSERVATION_DATA, PATH_CALCULATION_DATA, PATH_RESULT
from utilities import calculate_initial_coordinate_velocity
from date_helper import DateHelper
from file_helper import FileHelper

logger = Logger(PATH_ORBIT_LOG)
dh = DateHelper()
fh = FileHelper()

class OrbitalImprovement:
    def __init__(self):
        # Constants paremeter
        self.variation_coordinate = 10e-5
        self.variation_velocity = 10e-5
        self.variation_j20 = 10e-8
        self.variations = np.array([
            [self.variation_coordinate, 0e0, 0e0, 0e0, 0e0, 0e0],
            [0e0, self.variation_coordinate, 0e0, 0e0, 0e0, 0e0],
            [0e0, 0e0, self.variation_coordinate, 0e0, 0e0, 0e0],
            [0e0, 0e0, 0e0, self.variation_velocity, 0e0, 0e0],
            [0e0, 0e0, 0e0, 0e0, self.variation_velocity, 0e0],
            [0e0, 0e0, 0e0, 0e0, 0e0, self.variation_velocity],
            [0e0, 0e0, 0e0, 0e0, 0e0, 0e0]
        ])
        self.t = T / N
        self.initial_coordinate, self.initial_velocity = self.calculate_initial_coordinate_velocity()
        self.calculate_observation_data()
        self.delta_coordinate, self.delta_velocity, self.delta_j20 = self.generate_error()
        self.intermediate_coordinate = self.initial_coordinate
        self.intermediate_velocity = self.initial_velocity
        self.intermediate_j20 = J20
        
    def generate_error(self):
        # Add error to the initial coordinates and velocities
        modulus_coordinate = np.linalg.norm(self.initial_coordinate)
        modulus_velocity = np.linalg.norm(self.initial_velocity)
        
        arcsecond_radian = ARCSECOND * np.pi / 648000
        
        delta_coordinate = np.array([
            np.random.uniform(-1, 1) * arcsecond_radian * modulus_coordinate,
            np.random.uniform(-1, 1) * arcsecond_radian * modulus_coordinate,
            np.random.uniform(-1, 1) * arcsecond_radian * modulus_coordinate
        ])
        delta_velocity = np.array([
            np.random.uniform(-1, 1) * arcsecond_radian * modulus_velocity,
            np.random.uniform(-1, 1) * arcsecond_radian * modulus_velocity,
            np.random.uniform(-1, 1) * arcsecond_radian * modulus_velocity
        ])
        delta_j20 = np.random.uniform(-1, 1) / 100 * J20
        logger.log_info(f"Delta coordinate: {delta_coordinate}.")
        logger.log_info(f"Delta velocity: {delta_velocity}.")
        
        return delta_coordinate, delta_velocity, delta_j20

    # Get the initial coordinate and velocity
    def calculate_initial_coordinate_velocity(self):
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
        logger.log_info(f"Finished calculating intial coordinates and velocities. Elapsed time: {runtime:.6f} s.")
        logger.log_info(f"Initial coordinates: {initial_coordinate}, initial velocities: {initial_velocity}")
        print(f"Finished calculating intial coordinates and velocities. Elapsed time: {runtime:.6f} s.")
        print(f"Initial coordinates: {initial_coordinate}, initial velocities: {initial_velocity}")
        return initial_coordinate, initial_velocity

    # Compute the observation data using the integrator
    def calculate_observation_data(self):
        logger.log_info("Calculating observed coordinates and velocities.")
        print("Calculating observed coordinates and velocities. Please wait.")

        # Setup dataframe for saving data
        df = pd.DataFrame({
            "T0": [T0],
            "t": [0],
            "x_1": [self.initial_coordinate[0]],
            "x_2": [self.initial_coordinate[1]],
            "x_3": [self.initial_coordinate[2]],
            "vx_1": [self.initial_velocity[0]],
            "vx_2": [self.initial_velocity[1]],
            "vx_3": [self.initial_velocity[2]]
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
        # Set-up data
        t_i = dh.add_integration_time(self.t * 100)
        date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_i)

        # Configure and run the integrator
        fh.configure_file(
            coordinate=self.initial_coordinate,
            velocity=self.initial_velocity,
            date=date, 
            month=month, 
            year=year, 
            hour=hour, 
            minute=minute, 
            second=second, 
            t=self.t
        )
        fh.run_exe_file()
        
        for i in range(1, 101):
            logger.log_info(f"Finish calculating for observation data iteration: {i}.")
            logger.log_info(f"T0: {T0}")
            logger.log_info(f"t: {self.t * i}")
            logger.log_info(f"JD: ( {year} {month} {date} {hour} {minute} {second} )")
            logger.log_info(f"Computed coordinate: {fh.get_coordinate(i)}")
            logger.log_info(f"Computed velocity: {fh.get_velocity(i)}")

        # Get computed coordinates and velocities
        computed_coordinates, computed_velocities = fh.get_coordinates(), fh.get_velocities()

        # Add the computed coordinates and velocities to df
        computed_coordinates_df = pd.DataFrame(computed_coordinates, columns=["x_1", "x_2", "x_3"])
        computed_velocities_df = pd.DataFrame(computed_velocities, columns=["vx_1", "vx_2", "vx_3"])
        computed_df = pd.DataFrame(
            {
                "T0": [T0 for _ in range(1, 101)],
                "t": [self.t*i for i in range(1, 101)]
            }
        )
        computed_df = pd.concat([computed_df, computed_coordinates_df, computed_velocities_df], axis=1)
        df = pd.concat([df, computed_df], ignore_index=True)

        # Save the value to csv "observation_data.csv"
        df.to_csv(PATH_OBSERVATION_DATA, index=False)
        logger.log_info(f"Data is saved to {PATH_OBSERVATION_DATA}\n\n")
        print(f"Data is saved to {PATH_OBSERVATION_DATA}!")
        end = time()
        runtime = end - start
        logger.log_info(f"Finished calculating integration for intial coordinates and velocities. Elapsed time: {runtime:.6f} s.")
        print(f"Finished calculating integration for intial coordinates and velocities. Elapsed time: {runtime:.6f} s.")

    # Compute the calculated data using the integrator
    def calculate_calculated_data(self, n_iter):
        logger.log_info("Calculating computed coordinates and velocities.")
        print("Calculating computed coordinates and velocities. Please wait.")
        # Setup dataframe for saving data
        df = pd.DataFrame({
            "T0": [T0],
            "t": [0],
            "x_1": [self.initial_coordinate[0]],
            "x_2": [self.initial_coordinate[1]],
            "x_3": [self.initial_coordinate[2]],
            "vx_1": [self.initial_velocity[0]],
            "vx_2": [self.initial_velocity[1]],
            "vx_3": [self.initial_velocity[2]]
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
        
        error_coordinate = self.intermediate_coordinate + self.delta_coordinate
        error_velocity = self.intermediate_velocity + self.delta_velocity
        error_J20 = self.intermediate_j20 + self.delta_j20
        self.intermediate_coordinate = error_coordinate
        self.intermediate_velocity = error_velocity
        self.intermediate_j20 = error_J20
        logger.log_info(f"Error J20: {self.delta_j20}.")
        logger.log_info(f"Error coordinate: {self.delta_coordinate}.")
        logger.log_info(f"Error velocity: {self.delta_velocity}.")
        logger.log_info(f"J20 with error: {error_J20}.")
        logger.log_info(f"Coordinate with error: {error_coordinate}.")
        logger.log_info(f"Velocity with error: {error_velocity}.")

        # Generate the initial prediction for orbit (x10, x20, x30, x'10, x'20, x'30)
        # Set-up data
        t_i = dh.add_integration_time(self.t * 100)
        date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_i)

        # Configure and run the integrator
        fh.modify_garm360_in(new_value=error_J20)
        fh.configure_file(
            coordinate=error_coordinate,
            velocity=error_velocity,
            date=date, 
            month=month, 
            year=year, 
            hour=hour, 
            minute=minute, 
            second=second, 
            t=self.t
        )
        fh.run_exe_file()
        
        # for i in range(1, 101):
        #     logger.log_info(f"Finish calculating for calculated data iteration: {i}.")
        #     logger.log_info(f"T0: {T0}")
        #     logger.log_info(f"t: {self.t * i}")
        #     logger.log_info(f"JD: ( {year} {month} {date} {hour} {minute} {second} )")
        #     logger.log_info(f"Computed coordinate: {fh.get_coordinate(i)}")
        #     logger.log_info(f"Computed velocity: {fh.get_velocity(i)}")

        # Get computed coordinates and velocities
        computed_coordinates, computed_velocities = fh.get_coordinates(), fh.get_velocities()

        # Add the computed coordinates and velocities to df
        computed_coordinates_df = pd.DataFrame(computed_coordinates, columns=["x_1", "x_2", "x_3"])
        computed_velocities_df = pd.DataFrame(computed_velocities, columns=["vx_1", "vx_2", "vx_3"])
        computed_df = pd.DataFrame(
            {
                "T0": [T0 for _ in range(1, 101)],
                "t": [self.t*i for i in range(1, 101)]
            }
        )
        computed_df = pd.concat([computed_df, computed_coordinates_df, computed_velocities_df], axis=1)
        df = pd.concat([df, computed_df], ignore_index=True)

        # Save the value to csv "calculation_data.csv"
        df.to_csv(PATH_CALCULATION_DATA.format(n_iter), index=False)
        logger.log_info(f"Data is saved to {PATH_CALCULATION_DATA.format(n_iter)}\n\n")
        print(f"Data is saved to {PATH_CALCULATION_DATA.format(n_iter)}!")
        end = time()
        runtime = end - start
        logger.log_info(f"Finished calculating integration for intial coordinates and velocities. Elapsed time: {runtime:.6f} s.")
        print(f"Finished calculating integration for intial coordinates and velocities. Elapsed time: {runtime:.6f} s.")

    # Calcualte c vector
    def calculate_c(self, n_iter):
        print("Calculating difference between observation and computed data. Please wait.")
        logger.log_info("Calculating difference between observation and computed data. Please wait.")
        start = time()
        df_observation = pd.read_csv(PATH_OBSERVATION_DATA)
        df_calculation = pd.read_csv(PATH_CALCULATION_DATA.format(n_iter))

        c_df = df_observation[['x_1', 'x_2', 'x_3']] - df_calculation[['x_1', 'x_2', 'x_3']]
        c_vector = c_df.iloc[1:].values.flatten()
        end = time()
        runtime = end - start
        print(f"Finish calculating c vector. Elapsed time: {runtime:.6f} s.")
        logger.log_info(f"Finish calculating c vector: {c_vector}. Elapsed time: {runtime:.6f} s.")
        return c_vector

    # Calculate A (Isochronous derivatives matrix)
    def calculate_A(self, n_iter):
        print("Calculating A (Isochronous derivatives matrix). Please wait.")
        logger.log_info("Calculating A (Isochronous derivatives matrix).")
        start = time()
        df_calculation = pd.read_csv(PATH_CALCULATION_DATA.format(n_iter))

        coordinate_columns = ['x_1', 'x_2', 'x_3']

        A_matrix = []

        for column in range(7):
            A_column = []
            start_col = time()
            t_n = dh.add_integration_time(self.t * 100)
            date, month, year, hour, minute, second = dh.julian_date_to_gregorian(t_n)

            coordinate_velocity_with_variation = np.concatenate([
                self.initial_coordinate, 
                self.initial_velocity
            ]) + self.variations[column]
            coordinate_with_variation = np.split(coordinate_velocity_with_variation, 2)[0]
            velocity_with_variation = np.split(coordinate_velocity_with_variation, 2)[1]
            
            if column == 6:
                fh.modify_garm360_in(new_value=self.variation_j20)
            
            fh.configure_file(
                coordinate=coordinate_with_variation,
                velocity=velocity_with_variation,
                date=date, 
                month=month, 
                year=year, 
                hour=hour, 
                minute=minute, 
                second=second, 
                t=self.t
            )
            fh.run_exe_file()
            
            for i in range(1, 101):               
                for component in range(3):
                    computed_coordinate = df_calculation.iloc[i][coordinate_columns[component]]
                    computed_coordinate_variation = fh.get_coordinate(i)
                    
                    if column < 3:
                        variation = self.variation_coordinate
                    elif column >= 3 and column < 6:
                        variation = self.variation_velocity
                    else:
                        variation = self.variation_j20
                    derivative = (computed_coordinate_variation[component] - computed_coordinate) / variation
                    
                    A_column.append(derivative)
                
                logger.log_info(f"Calculating derivatives for row {i} at column {column + 1}: {derivative}.")
            A_matrix.append(A_column)
            end_col = time()
            runtime_col = end_col - start_col
            print(f"Finish calculating A matrix for column {column + 1}. Elapsed time: {runtime_col:.6f} s.")
            logger.log_info(f"Finish calculating A matrix for column: {column + 1}. Elapsed time: {runtime_col:.6f} s.")

        A_matrix = np.array(A_matrix).T

        end = time()
        runtime = end - start
        print(f"Finish calculating A Matrix. Elapsed time: {runtime:.6f} s.")
        logger.log_info(f"Finish calculating A Matrix with shape: {A_matrix.shape}. Elapsed time: {runtime:.6f} s.\n{A_matrix}")
        return A_matrix

    # Calculate Q
    def calculate_Q(self, A):
        print("Calculating Q. Please wait.")
        logger.log_info("Calculating y (desired correction).")
        Q = A.T @ A
        logger.log_info(f"Calculated Q with shape {Q.shape}: {Q}.")
        return Q

    # Calculate d
    def calculate_d(self, A, c):
        print("Calculating d. Please wait.")
        d = A.T @ c
        logger.log_info(f"Calculated d, with shape {d.shape}: {d}.")
        return d

    # Calculate y (desired correction)
    def calculate_y(self, Q, d):
        print("Calculating y (desired correction). Please wait.")
        y = np.linalg.inv(Q) @ d
        print(f"Finish calculating calculating y (desired correction).")
        logger.log_info(f"Finish calculating y (desired correction): {y}.")
        return y

    # Calculate sigma
    def calculate_sigma(self, Q, c):
        start = time()
        # Calculate root mean square error weight
        # Sigma0 = sqrt(a/b)
        # a = sum((delta x_ik(O-C))^2)
        # b = 3N-6
        logger.log_info(f"Calculating sigma 0 and sigma i.")
        print(f"Calculating sigma 0 and sigma i. Please wait.")
        a = np.sum(np.square(c))
        b = 3 * N - 7
        sigma_0 = np.sqrt(a/b)
        logger.log_info(f"Finish calculating Sigma 0: {sigma_0}")
        print(f"Sigma 0: {sigma_0}")

        # Calculate root mean square error parameters
        # Sigma i = Sigma 0 q_ii
        # q_ii = Diagonal element of Q^-1
        sigma_i = sigma_0 * np.diagonal(np.linalg.inv(Q))
        end = time()
        runtime = end - start
        
        logger.log_info(f"Sigma i: {sigma_i}")
        logger.log_info(f"Finish calculating Sigma i: {sigma_i}, norm: {np.linalg.norm(sigma_i)}. Elapsed time: {runtime:.6f} s.")
        print(f"Finished calculating sigma i: {sigma_i}, norm: {np.linalg.norm(sigma_i)}.")
        return sigma_0, sigma_i
    
    def run(self):
        # Prepare DataFrame to save the result
        iteration_label = "Iteration"
        sigma_label = "Sigma 0"
        sigma_x_label = "Sigma x"
        sigma_y_label = "Sigma y"
        sigma_z_label = "Sigma z"
        sigma_vx_label = "Sigma vx"
        sigma_vy_label = "Sigma vy"
        sigma_vz_label = "Sigma vz"
        sigma_j20_label = "Sigma J20"
        x_label = "X (km)"
        y_label = "Y (km)"
        z_label = "Z (km)"
        vx_label = "Vx (km/s)"
        vy_label = "Vy (km/s)"
        vz_label = "Vz (km/s)"
        J20_label = "J20"
        df = pd.DataFrame({
            iteration_label: [],
            sigma_label: [],
            J20_label: [],
            x_label: [],
            y_label: [],
            z_label: [],
            vx_label: [],
            vy_label: [],
            vz_label: [],
            sigma_x_label: [],
            sigma_y_label: [],
            sigma_z_label: [],
            sigma_vx_label: [],
            sigma_vy_label: [],
            sigma_vz_label: [],
            sigma_j20_label: []
        })
        
        print("Running main algorithm.")
        logger.log_info("Running main algorithm.")
        start = time()
        n_iter = 0
        prev_sigma_0 = float('inf')  # Start with an infinitely large sigma_0

        while True:
            n_iter += 1
            print(f"Iteration: {n_iter}.")
            logger.log_info(f"Iteration: {n_iter}.")

            # Compute calculation data
            self.calculate_calculated_data(n_iter)
            
            # Retrieve current coordinate and velocity
            coordinate = fh.get_coordinate(0)
            velocity = fh.get_velocity(0)
            
            # Calculate difference between observation and calculated, delta x_ik(O-C)
            c_vector = self.calculate_c(n_iter)
            
            # Calculate A (Isochronous derivatives matrix)
            A_matrix = self.calculate_A(n_iter)
            
            # Calculate y - desired correction
            Q = self.calculate_Q(A_matrix)
            d = self.calculate_d(A_matrix, c_vector)
            y = self.calculate_y(Q, d)
            
            self.delta_coordinate = y[:3]
            self.delta_velocity = y[3:6]
            self.delta_j20 = y[6]

            # Compute new sigma values
            sigma_0, sigma_i = self.calculate_sigma(Q, c_vector)

            # Store iteration data
            df.loc[len(df)] = [
                n_iter,
                sigma_0,
                self.intermediate_j20,
                coordinate[0],
                coordinate[1],
                coordinate[2],
                velocity[0],
                velocity[1],
                velocity[2],
                sigma_i[0],
                sigma_i[1],
                sigma_i[2],
                sigma_i[3],
                sigma_i[4],
                sigma_i[5],
                sigma_i[6]
            ]

            # Check if sigma_0 is still decreasing
            if sigma_0 >= prev_sigma_0:
                print(f"Stopping at iteration {n_iter} as sigma_0 is no longer decreasing.")
                break
            
            # Update previous sigma_0 for next iteration
            prev_sigma_0 = sigma_0
        
        end = time()
        runtime = end - start
        print(f"Finished running main algorithm. Elapsed time: {runtime:.6f} s.")
        logger.log_info(f"Finished running main algorithm. Elapsed time: {runtime:.6f} s.")
        df.to_csv(PATH_RESULT)
        
        
    
if __name__ == '__main__':
    om = OrbitalImprovement()
    om.run()