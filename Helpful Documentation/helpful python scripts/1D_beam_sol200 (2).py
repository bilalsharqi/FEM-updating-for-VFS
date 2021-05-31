"""
Quick and dirty script to produce beam SOL200 design bulk data cards (requires Python > 3.6 due to f-strings).
NOTE: This file requires user to delete a trailing + and duplicate row in objective function definition
"""

import math

# Give number of digits in number... will be used in writing objective function
def orderOfMag(a_number):
    return math.floor(math.log(a_number,10))+1


# -----------------------------------------------------------------
# ------------------------ File Definition ------------------------
# -----------------------------------------------------------------

manual_input = True
desvar_num = range(1,21) # Twenty design variables
height_ids = range(1, 21, 2)
width_ids = range(2,22,2)
var_label = range(1,11,1)

# ------------------------ Define Constraints ------------------------
node_constrain = range(1,22,1)
constraint = [0, # Grid 1
              -5.476546E-05, # Grid 2
              -2.092557E-04, # Grid 3
              -4.532520E-04, # Grid 4
              -7.770882E-04, # Grid 5
              -1.171650E-03, # Grid 6
              -1.628377E-03, # Grid 7
              -2.139259E-03, # Grid 8
              -2.696839E-03, # Grid 9
              -3.294213E-03, # Grid 10
              -3.925029E-03, # Grid 11
              -4.583487E-03, # Grid 12
              -5.264339E-03, # Grid 13
              -5.962891E-03, # Grid 14
              -6.675000E-03, # Grid 15
              -7.397076E-03, # Grid 16
              -8.126080E-03, # Grid 17
              -8.859527E-03, # Grid 18
              -9.595484E-03, # Grid 19
              -1.033257E-02, # Grid 20
              -1.106996E-02] # Grid 21

tolerance = 0.01 # +/- true value

# ------------------------ Define Design Space ------------------------
height_true = 0.0031 # Real value is 0.0031
height_lb = height_true*(1-tolerance)
height_ub = height_true*(1+tolerance)

width_true = 0.051 # Real value is 0.051
width_lb = width_true*(1-tolerance)
width_ub = width_true*(1+tolerance)

rho_init = 2810.0
rho_ub = 2810.0
rho_lb = 2640.0

if(manual_input):
    # Initial heights and width manually entered
    height_init = [height_lb,height_lb,height_lb,height_lb,height_lb,height_lb,height_lb,height_lb,height_lb,height_lb]
    width_init = [width_lb,width_lb,width_lb,width_lb,width_lb,width_lb,width_lb,width_lb,width_lb,width_lb]
else:
    print('Manual input false!')
    # Initial heights and widths random generated from range
    # ------



# ------------------------------------------------------------------
# ------------------------ Write Input File ------------------------
# ------------------------------------------------------------------
with open("sol200_problem.bdf", "w") as input_file:

# ------------------------ Write Height Design Variable ------------------------
    input_file.write('$...Height Design Variables\n')
    input_file.write('$ DESVAR, ID, LABEL, XINIT, XLB, XUB, DELXV,\n')
    for (des_id, label_num, init_val) in zip (height_ids, var_label, height_init):
        input_file.write(f'DESVAR,{des_id:8d},H_{label_num:d},{init_val:.2e},{height_lb:.2e},{height_ub:.2e},\n')
    input_file.write("\n")

# ------------------------ Write Width Design Variable ------------------------
    input_file.write('$...Width Design Variables\n')
    for (des_id, label_num, init_val) in zip (width_ids, var_label, width_init):
        input_file.write(f'DESVAR,{des_id:8d},W_{label_num:d},{init_val:.2e},{width_lb:.2e},{width_ub:.2e},\n')
    input_file.write("\n")
    
# ------------------------ Link Design Variable to PBEAML ------------------------
    input_file.write('$...Link Designable Properties to Design Variables. Grouped by BEAML\n')
    input_file.write('$ DVPREL1, ID, TYPE, PID, PNAME, PMIN, PMAX, C0,, ,\n')
    input_file.write('$        , DVID1, COEF1,\n')
    for (height_id,var,width_id) in zip(height_ids,var_label,width_ids):
        input_file.write(f'DVPREL1,{height_id*100},PBEAML,{var}, DIM2(A),\n')
        input_file.write(f'       ,{height_id:8d}, 1.0,\n')
        input_file.write(f'DVPREL1,{width_id*100},PBEAML,{var}, DIM1(A),\n')
        input_file.write(f'       ,{width_id:8d}, 1.0,\n')

# ------------------------ Write Design Response, Constraint, ObjFunction ------------------------
    input_file.write("\n$ Design response is out-of-plane deflection of select beam nodes\n")
    for(label_id, node) in zip(node_constrain,node_constrain):
        input_file.write(f'DRESP1,{label_id+150},D{node},DISP,,,3,,{node}\n')
    input_file.write('\n')

    # Lower and upper bounds flipped if constraint is positive
    input_file.write("$ Response constrained to be within tolerance of true value\n")
    for(true_value,node_num) in zip(constraint,node_constrain):
        input_file.write(f'DCONSTR,{node_num+250},{node_num+150},{true_value*(1+tolerance):.2e},{true_value*(1-tolerance):.2e}\n')
    input_file.write('\n')

    input_file.write("DCONADD,2000,")
    for(node_num,count) in zip(node_constrain,range(1,len(node_constrain)+1)):
        input_file.write(f'{node_num+250},')
        if(count % 7 == 0):
            input_file.write('+\n+,,')
    input_file.write('\n\n')

    # Write first line of DEQATN card 
    input_file.write("$ Objective Function pushes desvar to true value\n")
    input_file.write("$ ------------------ IMPORTANT: Delete duplicate row F1, trailing comma and + in last line ------------------ \n")
    input_file.write("DEQATN  21      F1(")
    chars_used = 4
    for(h_desvar,w_desvar) in zip(height_ids, width_ids):
        input_file.write(f'H{h_desvar},W{w_desvar},')
        chars_used += 6+orderOfMag(h_desvar)+orderOfMag(w_desvar)
        if(chars_used > (56 -(6+orderOfMag(h_desvar)+orderOfMag(w_desvar)) )):
            chars_used  = 0
            input_file.write("\n        ")
    input_file.write(')')
    
    # Assumes at least one design variable exists
    input_file.write("\n        ")
    input_file.write(f'= (H{height_ids[0]} - {height_true})**2;\n')

    for(h_desvar,w_desvar) in zip(height_ids,width_ids):
        input_file.write(f'        F{h_desvar} = (H{h_desvar} - {height_true})**2;\n')
        input_file.write(f'        F{w_desvar} = (W{w_desvar} - {width_true})**2;\n')
                         
                                  
        
    input_file.write("        R = ")
    chars_used = 4
    for(h_desvar,w_desvar) in zip(height_ids, width_ids):
        input_file.write(f'F{h_desvar}+F{w_desvar}+')
        chars_used += 6+orderOfMag(h_desvar)+orderOfMag(w_desvar)
        if(chars_used > (56 -(6+orderOfMag(h_desvar)+orderOfMag(w_desvar)) )):
            chars_used  = 0
            input_file.write("\n        ")
    input_file.write('\n\n')

# ------------------------ Write Final Design Response ------------------------
    input_file.write('$ Design response connects DESOBJ to DESVAR\n')
    input_file.write('DRESP2,22,R_1,21,,,,,,+\n')
    input_file.write('+,       DESVAR,')
    for(id,count) in zip(desvar_num,range(1,len(desvar_num)+1)):
        input_file.write(f'{id},')
        if(count % 7 == 0):
            input_file.write('+\n+,       ,,')

    input_file.write('\n\n$...Density Design Variables\n')
    for (des_id,mat_num) in zip(range(30,41),range(1,11)):
        input_file.write(f'DESVAR,{des_id:8d},D{des_id},{rho_init},{rho_lb},{rho_ub},\n')
        input_file.write(f'DVMREL1,{des_id+50},MAT1,{mat_num},RHO,\n')
        input_file.write(f'       ,{des_id},1.0,\n')
    input_file.write("\n")

    
        


    
        
        
