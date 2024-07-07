import spiceypy as spice
from datetime import datetime, timedelta
import numpy as np
import requests
import re
import sys

program_version = "2.2.0"
authors = ["H. A. Guler"]
verbose = False
AU = 1.495979e8 # km

help_keys = ["-h", "--help", "help"]
license_keys = ["--license", "license"]
citation_keys = ["--citation", "citation"]

help_text = """
MINORBIT2 Help
======================
MINORBIT2 is an n-body propagator for minor planets.

The inputs can be defined using a TXT file. Create a file containing the following
lines:

T0 YYYY-MM-DD
TF YYYY-MM-DD
DT XX

The T0 line defines the simulation initial time, TF line defines the simulation
end time, DT line defines the orbit propagator time step length in days. Then,
define the minor planets that will be simulated, with one minor planet on each line,
using MP lines as follows.

MP 2017 BX232
MP 2017 AC64
MP 2017 BM230

Save the TXT file and provide it as an argument to the program, i.e.:
'minorbit2.py example.txt'

If no files are provided as an argument, the program will ask for a filename at
initialization. If still no filename is provided (Enter is pressed), it will,
by default, attempt to read 'minorbit.txt'.

COPYRIGHT AND CITATION
======================
Minorbit2 is licensed under the terms of MIT License.
Please run the program with '--license' argument to print more details.

If you use MINORBIT2 or any derived software for your research, I'd like to
hear about it!

Please run the program with '--citation' argument to print a citation example.

ISSUES
======================
Bug reports should go to github.com/arda-guler/minorbit2
"""

citation_text = """
Example APA citation:

Guler, H. A., (2024). Minorbit2 (Version 2.2.0) [Source code].
Retrieved from https://github.com/arda-guler/minorbit2
"""

license_text = """
Minorbit2 is licensed under the terms of MIT License.
See the license text at https://github.com/arda-guler/minorbit2/blob/master/LICENSE
"""

class MainBody:
    def __init__(self, name, pos, GM):
        self.name = name
        self.pos = pos
        self.GM = GM # a.k.a. mu, standard gravitational param.

class MP:
    def __init__(self, des, pos, vel):
        self.des = des
        self.pos = pos
        self.vel = vel

def vprint(*i):
    global verbose
    if verbose:
        print(i)

def parseMPJPL(input_string):
    position_pattern = re.compile(r"X\s*=\s*([-+]?[\d.E+-]+)\s*Y\s*=\s*([-+]?[\d.E+-]+)\s*Z\s*=\s*([-+]?[\d.E+-]+)")
    velocity_pattern = re.compile(r"VX\s*=\s*([-+]?[\d.E+-]+)\s*VY\s*=\s*([-+]?[\d.E+-]+)\s*VZ\s*=\s*([-+]?[\d.E+-]+)")

    position_matches = position_pattern.findall(input_string)
    velocity_matches = velocity_pattern.findall(input_string)

    if len(position_matches) >= 2 and len(velocity_matches) >= 2:
        position_values = [float(value) for value in position_matches[1]]
        velocity_values = [float(value) for value in velocity_matches[1]]

        position_array = np.array(position_values)
        velocity_array = np.array(velocity_values)

        return [position_array, velocity_array]
    else:
        raise ValueError("The required lines were not found in the text")

def getMPJPL(designations, start_time):
    dt_stop_time = start_time + timedelta(days=1)
    stop_time = dt_stop_time.strftime("%Y-%m-%d")
    start_time = start_time.strftime("%Y-%m-%d")
    
    state_vectors = []
    
    for i, des in enumerate(designations):
        bodyid = "DES=" + des

        url = "https://ssd.jpl.nasa.gov/api/horizons.api"

        url += "?format=text&EPHEM_TYPE=VECTORS&OBJ_DATA=NO&VEC_TABLE=2&CENTER='500@0'"
        url += "&COMMAND='{}'&START_TIME='{}'&STOP_TIME='{}'&STEP_SIZE='3d'".format(bodyid, start_time, stop_time)

        response = requests.get(url)

        if response.status_code == 200:
            state_vectors.append(parseMPJPL(response.text))
            vprint(response.text)

        elif response.status_code == 400:
            print("Error 400 Bad Request")

        elif response.status_code == 405:
            print("Error 400 Method Not Allowed")

        elif response.status_code == 500:
            print("Error 400 Internal Server Error")

        elif response.status_code == 503:
            print("Error 400 Server Unavailable")

        if i > 1 and i % 20 == 0:
            print("Getting minor planet state vectors:", i / len(designations) * 100, "%")

    vprint("OK! Got minor planet state vectors from JPL!")
    return state_vectors

def grav_accel(mp, bodies):
    accel = np.array([0, 0, 0])
    
    for body in bodies:
        dist = np.linalg.norm(body.pos - mp.pos)
        grav_dir = (body.pos - mp.pos) / dist
        grav_mag = body.GM / dist**2
        
        accel = accel + grav_mag * grav_dir

    return accel

def stepYoshida8(mps, bodies, dt):
    # - - - CONSTANTS - - -
    w1 = 0.311790812418427e0
    w2 = -0.155946803821447e1
    w3 = -0.167896928259640e1
    w4 = 0.166335809963315e1
    w5 = -0.106458714789183e1
    w6 = 0.136934946416871e1
    w7 = 0.629030650210433e0
    w0 = 1.65899088454396

    ds = [w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7]

    cs = [0.3145153251052165, 0.9991900571895715, 0.15238115813844, 0.29938547587066, -0.007805591481624963,
          -1.619218660405435, -0.6238386128980216, 0.9853908484811935, 0.9853908484811935, -0.6238386128980216,
          -1.619218660405435, -0.007805591481624963, 0.29938547587066, 0.15238115813844, 0.9991900571895715,
          0.3145153251052165]
    # - - -   - - -   - - -

    for i in range(15):
        for mp in mps:
            mp.pos = mp.pos + mp.vel * cs[i] * dt

        accels = []
        for mp in mps:
            accels.append(grav_accel(mp, bodies))

        for imp, mp in enumerate(mps):
            mp.vel = mp.vel + accels[imp] * ds[i] * dt

    for mp in mps:
        mp.pos = mp.pos + mp.vel * cs[15] * dt

def read_config(sys_args):
    minute = 60
    hour = 60 * 60
    day = 60 * 60 * 24

    enter_to_quit = False
    
    if len(sys_args) < 2:
        enter_to_quit = True
        filename = input("Input filename (type 'help' for help): ")

        if not filename:
            filename = "minorbit.txt"

        else:
            if filename.lower() in help_keys:
                print(help_text)
                input("Press Enter to quit...")
                exit()

            elif filename.lower() in citation_keys:
                print(citation_text)
                input("Press Enter to quit...")
                exit()

            elif filename.lower() in license_keys:
                print(license_text)
                input("Press Enter to quit...")
                exit()

    else:
        filename = sys_args[1]

        if filename.lower() in help_keys:
            print(help_text)
            input("Press Enter to quit...")
            exit()

        elif filename.lower() in citation_keys:
            print(citation_text)
            input("Press Enter to quit...")
            exit()

        elif filename.lower() in license_keys:
            print(license_text)
            input("Press Enter to quit...")
            exit()
        
        if not "." in filename:
            filename = filename + ".txt"

    print("")
    print("Reading inputs from " + filename + "...")

    MPdes = []
    pattern = re.compile(r'(\d{4})\s([A-Za-z0-9]+)')

    while True:
        try:
            with open(filename, "r") as f:
                lines = f.readlines()

                for line in lines:
                    if line.startswith("T0 "):
                        t_0 = datetime.strptime(line.split(";")[0][3:13].replace(" ", "").replace(" ", ""), "%Y-%m-%d")
                    elif line.startswith("TF "):
                        t_f = datetime.strptime(line.split(";")[0][3:13].split(" ")[0].replace(" ", "").replace(" ", ""), "%Y-%m-%d")
                    elif line.startswith("DT "):
                        dt = float(line[3:].split(";")[0].split(" ")[0].replace(" ", "").replace(" ", "")) * day
                    elif line.startswith("MP "):
                        match = pattern.search(line)
                        MPdes.append(match.group(1) + " " + match.group(2))
                    elif line.startswith("RF "):
                        result_filename = line[3:].split(";")[0].strip()
                        
        except FileNotFoundError:
            print("File '" + filename + "' not found!")

            if not filename.endswith(".txt"):
                filename = filename + ".txt"
                continue
                
            else:
                if enter_to_quit:
                    input("Press Enter to quit.")

                quit()

        print("Loading inputs from", filename, "...")
        break

    if not MPdes:
        print("Warning: No minor planet designations found in input file.")

    print("")
    print("Inputs:")
    
    try:
        print("T0", t_0)
    except UnboundLocalError:
        print("ERROR: T0 not assigned!")
        if enter_to_quit:
            input("Press Enter to quit.")
        quit()
        
    try:
        print("TF", t_f)
    except UnboundLocalError:
        print("ERROR: TF not assigned!")
        if enter_to_quit:
            input("Press Enter to quit.")
        quit()

    try:
        print("DT", dt, "seconds")
    except UnboundLocalError:
        print("ERROR: DT not assigned!")
        if enter_to_quit:
            input("Press Enter to quit.")
        quit()
        
    try:
        print("RF", result_filename)
    except UnboundLocalError:
        print("ERROR: RF not assigned!")
        if enter_to_quit:
            input("Press Enter to quit.")
        quit()

    if MPdes and len(MPdes) < 20:
        print("DES", end=" ")

        for d in MPdes[:-1]:
            print(d, end=", ")
        print(MPdes[-1])

    elif MPdes:
        print("There are lots of designations. Not printing them all for performance.")

    elif len(MPdes) == 0:
        print("INFO: No minor planets found in input file.")

    print("")
    
    return t_0, t_f, dt, MPdes, result_filename, enter_to_quit

def main(sys_args):
    global verbose
    
    print(" = = = MINORBIT = = =")
    print("Version", program_version)
    print("Authors: ", end="")

    for a in authors[:-1]:
        print(a, end=", ")
    print(authors[-1])
        
    print("")
    
    ## = = = SETUP = = =
    t_0, t_f, dt, MP_des, result_filename, enter_to_quit = read_config(sys_args)
    date_init = t_0
    date_final = t_f

    spice.furnsh('data/naif0012.tls')
    spice.furnsh('data/de421.bsp')

    mps = []
    MP_statevectors = getMPJPL(MP_des, date_init)
    for i in range(len(MP_des)):
        sv = MP_statevectors[i]
        new_MP = MP(MP_des[i], sv[0], sv[1])
        mps.append(new_MP)

    time_interval = (date_final - date_init).total_seconds()
    N_cycles = int(time_interval // dt) + 1
    date_final_actual = date_init + timedelta(seconds=N_cycles * dt)

    bodies = []
    body_names = ["MERCURY BARYCENTER",
                  "VENUS BARYCENTER",
                  "EARTH BARYCENTER",
                  "MARS BARYCENTER",
                  "JUPITER BARYCENTER",
                  "SATURN BARYCENTER",
                  "URANUS BARYCENTER",
                  "NEPTUNE BARYCENTER",
                  "SUN"]
    
    body_GMs = [2.2031780000000021E+04,
                3.2485859200000006E+05,
                4.0350323550225981E+05,
                4.2828375214000022E+04,
                1.2671276480000021E+08,
                3.7940585200000003E+07,
                5.7945486000000080E+06,
                6.8365271005800236E+06,
                1.3271244004193938e11]

    for i in range(9):
        new_body = MainBody(body_names[i], np.array([0, 0, 0]), body_GMs[i])
        bodies.append(new_body)

    MP_poses = []
    result_str = ""

    # == == MAIN LOOP == ==
    print("")
    print("Propagation starts...")
    for cycle in range(N_cycles):
        # get cycle datetime
        cycle_date = date_init + timedelta(seconds=cycle * dt)
        cycle_date = cycle_date.strftime('%Y-%m-%dT%H:%M:%S')
        t = spice.str2et(cycle_date)
        
        # update main body positions
        for ib, body in enumerate(bodies):
            state, _ = spice.spkezr(body.name, t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
            body.pos = state[:3]

        # update minor planet state vectors
        stepYoshida8(mps, bodies, dt)

        # update results
        result_str += "\n" + str(cycle_date) + "\t" + str(cycle * dt) + "\t"
        for b in bodies:
            result_str += str(','.join(map(str, b.pos))) + "\t"

        for mp in mps:
            result_str += str(','.join(map(str, mp.pos))) + "\t"

        # print progress
        if cycle % 200 == 0:
            print(str(cycle / N_cycles * 100) + "%")

        # CUSTOM PLOTS AND DATA...
        # lorem.ipsum(dolor, sit, amet)
        # sari = cizmeli(mehmet, aga)
        # ...
        
    print("Propagation complete.")
    # == == END MAIN LOOP == ==

    # check accuracy
    mp_pos_errors = []
    MP_statevectors_final = getMPJPL(MP_des, date_final_actual)
    for i, mp in enumerate(mps):
        sv = MP_statevectors_final[i]
        sv_pos = sv[0]
        error = np.linalg.norm(mp.pos - sv[0]) / AU
        mp_pos_errors.append(error)

    print("")
    print("Validation - Minor planet final position errors:")
    for error in mp_pos_errors:
        print(error, "AU")

    # file output
    print("")
    print("Writing results to", result_filename, "...")
    with open(result_filename, "w") as f:
        f.write(result_str)

    print("Results exported.")
    
    spice.unload('naif0012.tls')
    spice.unload('de421.bsp')

    print("")
    print("END PROGRAM")

    if enter_to_quit:
        input("Press Enter to quit.")

if __name__ == "__main__":
    main(sys.argv)
