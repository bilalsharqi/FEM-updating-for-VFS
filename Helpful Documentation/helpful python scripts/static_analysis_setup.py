"""
Quick and dirty script to produce uCRM SOL200 design bulk data cards (requires Python > 3.6 due to f-strings).
NOTE: This file requires user to delete a trailing + and duplicate row in objective function definition
"""

import math
from random import seed
from random import random

# Give number of digits in number... will be used in writing objective function
def orderOfMag(a_number):
    return math.floor(math.log(a_number,10))+1


# -----------------------------------------------------------------
# ------------------------ File Definition ------------------------
# -----------------------------------------------------------------

manual_input = False
desvar1_ids = range(1, 101, 1)

# ------------------------ Define Constraints ------------------------

# Constraint is modal frequency
constraints = []

with open("reference_results.txt") as ref_result:
    for line in ref_result:
       constraints.append(float(line))
    
constraint_numbers = range(1,len(constraints)+1)

tolerance = 0.05 # +/- true value

# ------------------------ Define Design Space ------------------------

height_true = 0.031
height_lb = height_true*(1-tolerance)
height_ub = height_true*(1+tolerance)

width_true = 0.51
width_lb = width_true*(1-tolerance)
width_ub = width_true*(1+tolerance)

# ------------------------ Define Structure Grids ------------------------
beam_length = 15 # meters
beam_dis = 100

with open("mat_and_props_BEAM.bdf","w") as model_file:
    model_file.write("$ Nodes for the entire model.\n")
    for grid_id in range(0,beam_dis+1):
        model_file.write(f'GRID     {grid_id+1}              {(beam_length/beam_dis)*grid_id:.2f}      0.      0.\n')

    model_file.write('\n$ Material definition for each beam section\n')
    model_file.write('$ Assumes material is aluminum\n')
    for mat_id in range(1,beam_dis+1):
        model_file.write(f'MAT1     {mat_id}      7.+10           .33     2700.\n')
        
    model_file.write('\n$ PBEAML definition for each section of beam\n')
    for pbeam_id in range(1,beam_dis+1):
        model_file.write(f'PBEAML   {pbeam_id}       {pbeam_id}               BAR\n')
        model_file.write('        .51     .031\n')

    model_file.write('\n$ CBEAM definition along beam length\n')
    for cbeam_id in range(1,beam_dis+1):
        model_file.write(f'CBEAM,{cbeam_id},{cbeam_id},{cbeam_id},{cbeam_id+1},0.,1.,0.\n')


seed(1)
with open("sol200_problem.bdf","w") as input_file:
    
# ------------------------ Write Height Design Variable ------------------------
    input_file.write('$...Height Design Variables\n')
    input_file.write('$ DESVAR, ID, LABEL, XINIT, XLB, XUB, DELXV,\n')
    for des_id in range(1,beam_dis+1):
        random_height = height_lb + random()*(height_ub - height_lb)
        input_file.write(f'DESVAR,{des_id:8d},H_{des_id:d},{random_height:.2e},{height_lb:.2e},{height_ub:.2e},\n')
    input_file.write("\n")

# ------------------------ Write Width Design Variable ------------------------
    input_file.write('$...Width Design Variables\n')
    for des_id in  range(beam_dis+1,2*beam_dis+1):
        random_width = width_lb + random()*(width_ub - width_lb)
        input_file.write(f'DESVAR,{des_id:8d},W_{des_id:d},{random_width:.2e},{width_lb:.2e},{width_ub:.2e},\n')
    input_file.write("\n")
    
# ------------------------ Link Design Variable to PBEAML ------------------------
    input_file.write('$...Link Designable Properties to Design Variables. Grouped by BEAML\n')
    input_file.write('$ DVPREL1, ID, TYPE, PID, PNAME, PMIN, PMAX, C0,, ,\n')
    input_file.write('$        , DVID1, COEF1,\n')
    for (height_id,width_id) in zip(range(1,beam_dis+1),range(beam_dis+1,2*beam_dis+1)):
        input_file.write(f'DVPREL1,{height_id},PBEAML,{height_id}, DIM2(A),\n')
        input_file.write(f'       ,{height_id:8d}, 1.0,\n')
        input_file.write(f'DVPREL1,{width_id},PBEAML,{height_id}, DIM1(A),\n')
        input_file.write(f'       ,{width_id:8d}, 1.0,\n')
    
# ------------------------ Write Design Response, Constraint, ObjFunction ------------------------
    input_file.write("\n$ Design response is out-of-plane deflection of select beam nodes\n")
    for(label_id) in (constraint_numbers):
        input_file.write(f'DRESP1,{label_id},D{label_id},DISP,,,3,,{label_id}\n')
    input_file.write('\n')

    # Modal frequency is constrained here. Variable names can be adapted to any constraint type however.
    input_file.write("$ Response constrained to be within tolerance of true value\n")
    for(constraint_num, true_value) in zip(constraint_numbers, constraints):
        if(true_value > 0):
            input_file.write(f'DCONSTR,{constraint_num},{constraint_num},{true_value*(1-tolerance):.2e},{true_value*(1+tolerance):.2e}\n')
        elif(true_value < 0):
            input_file.write(f'DCONSTR,{constraint_num},{constraint_num},{true_value*(1+tolerance):.2e},{true_value*(1-tolerance):.2e}\n')
        else:
            input_file.write(f'DCONSTR,{constraint_num},{constraint_num},-0.001,0.001\n')     
    input_file.write('\n')

    input_file.write("DCONADD,2000,")
    for(constraint_num) in (constraint_numbers):
        input_file.write(f'{constraint_num},')
        if(constraint_num % 7 == 0):
            input_file.write('+\n+,,')
    input_file.write('\n\n')


 # Write first line of DEQATN card 
    input_file.write("$ Objective Function pushes desvar to true value\n")
    input_file.write("$ ------------------ IMPORTANT: Delete duplicate row (F_1) adnd trailing + in last line ------------------ \n")
    input_file.write("DEQATN  21      F1(")
    chars_used = 3
    for(h_desvar) in range(1,beam_dis+1):
        input_file.write(f'H{h_desvar},')
        chars_used += 2+orderOfMag(h_desvar)
        if(chars_used > (56 -(3+orderOfMag(h_desvar)) )):
            chars_used  = 0
            input_file.write("\n        ")
    input_file.write('\n')
    input_file.write("        W101,")
    chars_used = 4
    
    for(w_desvar) in range(beam_dis+2,2*beam_dis+1):
        input_file.write(f'W{w_desvar},')
        chars_used += 2+orderOfMag(w_desvar)
        if(chars_used > (56 -(3+orderOfMag(w_desvar)) )):
            chars_used  = 0
            input_file.write("\n        ")
    input_file.write(')')

    for(h_desvar) in range(1,beam_dis+1):
        input_file.write(f'        F{h_desvar} = 100.0*(H{h_desvar} - {height_true})**2;\n')
    for(w_desvar) in range(beam_dis+1,2*beam_dis+1):
        input_file.write(f'        F{w_desvar} = (W{w_desvar} - {width_true})**2;\n')
    
    # Return HERE
    input_file.write("        R = ")
    chars_used = 4
    for(desvar1_id) in range(1,2*beam_dis+1):
        input_file.write(f'F{desvar1_id}+')
        chars_used += 2+orderOfMag(desvar1_id)
        if(chars_used > (56 - (2+orderOfMag(desvar1_id)) )):
            chars_used  = 0
            input_file.write("\n        ")
    input_file.write('\n\n')

    # ------------------------ Write Final Design Response ------------------------
    input_file.write('$ Design response connects DESOBJ to DESVAR\n')
    input_file.write('DRESP2,1000,R_1,21,,,,,,+\n')
    input_file.write('+,       DESVAR,')
    for(desvar_id) in range(1,2*beam_dis+1):
        input_file.write(f'{desvar_id},')
        if(desvar_id % 7 == 0):
            input_file.write('+\n+,       ,,')
 
