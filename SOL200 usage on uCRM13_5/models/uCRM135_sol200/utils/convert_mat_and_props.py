"""
Quick and dirty script to modify uCRM materials/properties (requires Python > 3.6 due to f-strings).
"""


E = 72.0e9
nu = 0.33
rho = 2780.0

with open("../uCRM-135_mat_and_props_reference.bdf", "r") as ref_file:
    with open("../uCRM-135_mat_and_props.bdf", "w") as input_file:
        for line1 in ref_file:
            line2 = next(ref_file)
            line3 = next(ref_file)
            # Second MAT2*
            line4 = next(ref_file)
            line5 = next(ref_file)
            line6 = next(ref_file)
            # Third MAT2*
            line7 = next(ref_file)
            line8 = next(ref_file)
            line9 = next(ref_file)
            # Fourth MAT2*
            line10 = next(ref_file)
            line11 = next(ref_file)
            line12 = next(ref_file)
            # PSHELL*
            line13 = next(ref_file)
            line14 = next(ref_file)
            line15 = next(ref_file)
            # CORD2R
            line16 = next(ref_file)
            line17 = next(ref_file)

            # Get numbering and current thickness from PSHELL card.
            id = int(line13.split()[1])
            t = float(line13.split()[3])

            # Write new file.
            input_file.write(f'MAT1, {id}, {E:.2e}, , {nu:.2f}, {rho:.1f},\n')
            input_file.write(f'PSHELL, {id}, {id},{t:.2e}, {id}, , {id},\n')
            input_file.write(line16)
            input_file.write(line17)
