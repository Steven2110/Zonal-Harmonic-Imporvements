import numpy as np
import math

from constants import GM, INPUT_TEMPLATE

# Newton's Method
def solve_kepler_newton(M, e, tol=1e-15, max_iter=1000):
    E = M  # Initial guess
    for _ in range(max_iter):
        f_E = E - e * np.sin(E) - M  # f(E)
        f_prime_E = 1 - e * np.cos(E)  # f'(E)

        # Update E using the iterative formula
        E_next = E - f_E / f_prime_E
        # Check the convergence condition
        if abs(E - E_next) < tol:
            return E_next  # Accept solution if convergence criterion is met

        E = E_next  # Update E for the next iteration

    return E  # Return the final value if max_iter is reached

def calculate_initial_coordinate_velocity(
    semi_major_axis: float,
    eccentricity: float,
    longitude_of_ascending_node: float,
    argument_pericenter: float,
    inclination: float,
    mean_anomaly: float
):
    E = solve_kepler_newton(mean_anomaly, eccentricity)
    # Calculate the value of 𝛏 and 𝛈
    ksi = semi_major_axis * (np.cos(E) - eccentricity)
    eta = semi_major_axis * np.sqrt(1-eccentricity**2) * np.sin(E)

    # Calculate all sin & cos of some variables
    sin_omega = np.sin(argument_pericenter)
    cos_omega = np.cos(argument_pericenter)
    sin_Omega = np.sin(longitude_of_ascending_node)
    cos_Omega = np.cos(longitude_of_ascending_node)
    sin_i = np.sin(inclination)
    cos_i = np.cos(inclination)

    # Calculate coordinates position (x, y, z)
    x = ksi * (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i) + eta * (-sin_omega * cos_Omega - cos_omega * sin_Omega * cos_i)
    y = ksi * (cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i) + eta * (-sin_omega * sin_Omega + cos_omega * cos_Omega * cos_i)
    z = ksi * sin_omega * sin_i + eta * cos_omega * sin_i

    r_vec_inertial = np.array([x, y, z])

    # Calculate the p & r
    p = semi_major_axis * (1 - eccentricity**2)
    r = semi_major_axis * (1 - eccentricity * np.cos(E)) # or r = np.linalg.norm(r_vec_inertial)

    # Calculate u
    sin_u = z / (r * sin_i)
    cos_u = x / r * cos_Omega + y / r * sin_Omega

    # Calculate the true anomaly ν
    sin_nu = sin_u * cos_omega - cos_u * sin_omega
    cos_nu = cos_u * cos_omega + sin_u * sin_omega

    # Calculate Vr and Vn
    Vr = np.sqrt(GM/p) * eccentricity * sin_nu
    Vn = np.sqrt(GM/p) * (1 + eccentricity * cos_nu)

    # Calculate velocity
    v_x = x/r * Vr + (-sin_u * cos_Omega - cos_u * sin_Omega * cos_i) * Vn
    v_y = y/r * Vr + (-sin_u * sin_Omega + cos_u * cos_Omega * cos_i) * Vn
    v_z = z/r * Vr + cos_u * sin_i * Vn

    v_vec_inertial = np.array([v_x, v_y, v_z])

    return r_vec_inertial, v_vec_inertial

def format_scientific(num):
    if num == 0.0:
        return "0.00000000000000E+0000"  # Handle zero case separately

    exponent = math.floor(math.log10(abs(num)))  # Find exponent
    mantissa = num / (10**exponent)  # Normalize to 1 digit before decimal
    formatted_str = f"{mantissa:.14f}E{exponent:+05d}"  # Format output

    return formatted_str

def rewrite_input_file(q, hour, minute, second, t):
    path_input_file = "iszm_puc.in"
    f = open(path_input_file, "w", encoding="cp866")
    f.write(INPUT_TEMPLATE.format(
        format_scientific(q[0]),
        format_scientific(q[1]),
        format_scientific(q[2]),
        format_scientific(q[3]),
        format_scientific(q[4]),
        format_scientific(q[5]),
        hour,
        minute,
        second,
        t
    ))
    
def get_date(t):
    h = int(t // 3600)
    sec = t - h * 3600
    m = int(sec // 60)
    sec -= m * 60
    s = round(sec, 5)
    return h, m, s