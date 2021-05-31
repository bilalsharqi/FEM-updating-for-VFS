# MSC Nastran Interface


# Author: 
# Cristina Riso, criso@umich.edu


# python modules
 
import numpy as np
import os
import shutil


def runNastran(inp_dir, run_dir, out_dir, run_type, debug=True):
    """
    Runs a MSC Nastran analysis using pre-set input file
    
    Args:
        inp_dir: path to the directory with the input files
        run_dir: path to the directory where to run the analysis
        out_dir: path to the directory where to store the output files
        run_type: analysis type ("modal" or "static")
        debug: flag for debug messages
        
    Remarks:
        1. The MSC Nastran modal analysis MUST NOT ask for an .op2 output. 
           In that case, the while loop never terminates (needs cleaning). 
        2. The analysis driver musut be called sol103*.dat or sol101*.dat. 
        3. Output file names are determined accordingly.
        4. Input files must be available in the inp_dir directory.
    """
    
    fname = "pyMSCNastranUtilities.runNastran"
    
    if debug:
        print("Running MSC Nastran analysis...\n")
    
    # checking inputs
    
    if not os.path.isdir(inp_dir) or os.listdir(inp_dir) == []:        
        raise ValueError("Error in "+fname+": input files not found.")
        
    if not run_type in ["modal","static"]:
        raise ValueError("Error in "+fname+": invalid run type.")
    elif run_type == "modal":
        driver = "sol103.dat"
    elif run_type == "static":
        driver = "sol101.dat"
    
    if not os.path.isfile(os.path.join(inp_dir,driver)):
        raise ValueError("Error in "+fname+": analysis driver not found.")
    
    # creating run and output directories if not existing or cleaning
    
    for dir_name in [run_dir,out_dir]:
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        elif dir_name == run_dir:
            for file in os.listdir(run_dir):
                os.remove(os.path.join(run_dir,file))

    # copying input files in the run directory
    
    inp_files = os.listdir(inp_dir)

    for inp_file in inp_files:

        if inp_file == "__pycache__":
            continue

        old_path = os.path.join(inp_dir,inp_file)
        new_path = os.path.join(run_dir,inp_file)

        shutil.copy(old_path,new_path)

    # calling MSC Nastran and wait for the analysis to end

    old_cwd = os.getcwd()
    new_cwd = run_dir

    os.chdir(new_cwd)

    os.system('nastran '+driver)

    while os.path.isfile(driver.replace(".dat",".op2")):
        pass

    os.chdir(old_cwd)

    # moving inputs and outputs to the output directory
    
    out_files = []
    
    for file in os.listdir(run_dir):
        if ".MASTER" in file or ".DBALL" in file or ".IFPDAT" in file:
            os.remove(os.path.join(run_dir,file))
    
    for file in os.listdir(run_dir):
        if "sol" in file and not ".dat" in file:
            out_files.append(file)

    for file in inp_files+out_files:

        if file == "__pycache__":
            continue

        old_path = os.path.join(run_dir,file)
        new_path = os.path.join(out_dir,file)

        if os.path.isfile(new_path):
            os.remove(new_path)

        os.rename(old_path,new_path)
        
    for file in os.listdir(run_dir):
        os.remove(os.path(run_dir,file))

    if debug:
        print("\nMSC Nastran analysis completed\n")

    return


def importGrids(out_dir, out_files, debug=True):
    """
    Imports the FEM grid IDs and coordinates.
    
    Args:
        out_dir: path to the output directory        
        out_files: names of files containg the FEM grid IDs and coordinates
        debug: flag for debug messages        
        
    Returns:
        grid_ids: FEM grid IDs (list of integers)
        grid_coords = FEM grid coordinates (n_grids*3 numpy array)
        
    Remarks:
        1. Grid IDs are imported from out_files[0] which must be a BDF file.
        2. Grid coordinates are imported from out_files[1] which must be a DMIG
           format file obtained from, for instance, a SOL103 ALTER output.
        3. Grid IDs are returned in increasing order.
        4. Grid coordinates are returned in the global coordinate system.
    """
    
    fname = "pyMSCNastranUtilities.importGrids"
    
    # checking inputs
    
    if len(out_files) < 2:
        raise ValueError("Error in "+fname+": provide output file names.")

    # importing grid IDs and total number of model grids
    
    grid_ids, n_grids = importGridIDs(out_dir,out_files[0],debug)

    # importing grid coordinates

    grid_coords = importGridCoords(out_dir,out_files[1],grid_ids,debug)

    return grid_ids, n_grids, grid_coords


def importGridIDs(out_dir, out_file, debug=True):
    """
    Imports the FEM grid IDs from a BDF file.
    
    Args:
        out_dir: path to the output directory
        out_file: BDF file name
        debug: flag for debug messages        
        
    Returns:
        grids: grid IDs (list of integers)
        n_grids: number of grids (scalar)
        
    Remarks:
        1. Grid IDs must be in long/short or comma-separated formats.
        2. Grid IDs are returned in increasing order.
    """
    
    fname = "pyMSCNastranUtilities.importGridIDs"    

    if debug:
        print("Importing grids...")

    # checking inputs

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": BDF file not found.")

    # reading input file

    with open(os.path.join(out_dir,out_file),"r") as open_file:
        lines = open_file.readlines()
    open_file.close()
    
    # allocation
    
    grids = []

    # looking for GRID entries  

    for line in lines:

        if line.startswith("GRID"):

            # formatted short or long format

            if line[4] == " " or line[4] == '*':
                
                grid = line.split()[1]
                
                if line[7]=='*':
                
                    grid = line.split()[2]
                    
                if grid[0] == "*":

                    grid = grid[1:]

            # comma-separated format

            elif line[4] == ",":

                grid = line.split(",")[1]

            else:

                raise ValueError("Error in "+fname+": invalid GRID format.")

            # storing GRID ID
            
            grids.append(int(grid))
            
    if debug:
        print("Grids imported")            

    return grids, len(grids)


def importGridCoords(out_dir, out_file, grids, debug=True):
    """
    Imports the FEM grid coordinates from a DMIG format file. 
    
    Args:
        out_dir: path to the output directory
        out_file: DMIG format file name        
        grids: IDs of grids whose coordinates have to be imported
        debug: flag for debug messages
        
    Returns:
        x: grid coordinates (n_grids*3 numpy array)
        
    Remarks:
        1. The DMIG format file can be obtained from a SOL103 ALTER output.
        2. Grid are returned in increasing ID order.
        3. Grid coordinates are returned in the global coordinate system.
    """
    
    fname = "pyMSCNastranUtilities.importGridCoords"    

    if debug:
        print("Importing grids coordinates...")
        
    # checking inputs
    
    if grids == []:
        raise ValueError("Error in "+fname+": provide grid IDs.")

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")

    # reading file

    with open(os.path.join(out_dir,out_file),"r") as open_file:
        lines = open_file.readlines()
    open_file.close()
    
    # allocation

    x = np.zeros([3,len(grids)])    

    # importing grid coordinates

    for line in lines:

        if line.startswith("*"):

            line_split = line.split()

            grid = line_split[1]
            comp = line_split[2][0]

            if "-" in line_split[2]:

                value = line_split[2][1:].replace("D","e")

            else:

                value = line_split[3].replace("D","e")
                
            try:
                x[int(comp)-1,grids.index(int(grid))] = value
            except:
                raise ValueError("Error in "+fname+": grid "+grid+" not found.")
                
    if debug:
        print("Grid coordinates imported")  

    return x


def importLumpedMassData(out_dir, out_file, grids, debug=True):
    """
    Imports the FEM lumped mass data from a DMIG format file. 
    
    Args:
        out_dir: path to the output directory
        out_file: DMIG format file name        
        grid: FEM grid IDs whose lumped mass data has to be imported
        debug: flag for debug messages
        
    Returns:
        m: nodal masses (1*n_grids numpy array)
        o: nodal mass offsets (3*n_grids numpy array)
        J: nodal inertia tensors (list of 3*3 numpy arrays)
        
    Remarks:
        1. Data is returned in increasing grid ID order.
        2. Components are in the global coordinate system. 
    """
    
    fname = "pyMSCNastranUtilities.importLumpedMassData"

    if debug:
        print("Importing lumped mass data...")

    # checking inputs

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")

    # reading file

    with open(os.path.join(out_dir,out_file),"r") as open_file:
        lines = open_file.readlines()
    open_file.close()    

    # allocations
    
    n_grids = len(grids)

    m = np.zeros([1,n_grids])
    o = np.zeros([3,n_grids])
    J = np.zeros([6,n_grids])

    # importing lumped mass data

    for line in lines:

        if line.startswith("DMIG*"):

            line_split = line.split()

            grid = int(line_split[2])

            mcol = int(line_split[3])
            
            try:
                index = grids.index(grid)
            except:
                raise ValueError("Error in "+fname+": grid "+grid+" not found.")                

        elif line.startswith("*"):

            line_split = line.split()

            mrow = int(line_split[2][0])

            if "-" in line_split[2]:

                value = float(line_split[2][1:].replace("D","e"))

            else:

                value = float(line_split[3].replace("D","e"))            

            if mcol == 1:

                m_i = value

                m[0,index] = m_i

            elif mcol == 4:

                if mrow in [2, 3]:

                    o_i = value/m_i

                    if mrow == 2:

                        o[2,index] = -o_i

                    else:

                        o[1,index] = +o_i

                else:

                    J[0,index] = value

            elif mcol == 5:

                if mrow == 3:

                    o_i = value/m_i

                    o[0,index] = -o_i

                elif mrow in [4, 5]:

                    J[mrow-3,index] = value

            elif mcol == 6:

                if mrow in [4, 5, 6]:

                    J[mrow-1,index] = value

    # recasting J as list of inertia tensors

    J_new = list(range(0,n_grids))

    for i in range(0,n_grids):

        J_i = np.zeros([3,3])

        J_i[0,0] = J[0,i]
        J_i[0,1] = J[1,i]
        J_i[1,1] = J[2,i]
        J_i[0,2] = J[3,i]
        J_i[1,2] = J[4,i]
        J_i[2,2] = J[5,i]

        J_i[1,0] = J_i[0,1]
        J_i[2,0] = J_i[0,2]
        J_i[2,1] = J_i[1,2]

        J_new[i] = J_i
        
    if debug:
        print("Lumped mass data imported")

    return m, o, J_new


def importDisplacements(out_dir, out_file, n_subcases, grids, grids_order=[], debug=True):
    """
    Imports a set of displacement fields from a F06 file.
    
    Args:
        out_dir: path to the output directory
        out_file: F06 file name
        n_fields: number of displacement fields to be imported
        grids: IDs of the grids whose displacements have to be imported
        debug: flag for debug messages
        
    Returns:
        u: displacement fields (list of 6*n_grids numpy arrays)
        
    Remarks:
        1. Components are in the output coordinate system of each grid.
        2. Displacement fields are returned in provided grid ID order.
    """
    
    fname = "pyMSCNastranUtilities.importDisplacements"

    if debug:
        print("Importing static displacements...")
        
    # checking inputs
    
    if grids == []:
        raise ValueError("Error in "+fname+": provide grid IDs.")

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")        

    # allocations
    
    n_grids = len(grids)

    u = [np.zeros([6,n_grids]) for i in range(n_subcases)]

    # importing displacement fields

    open_file = open(os.path.join(out_dir,out_file),"r") 
    
    cnt = 0
    for count, line in enumerate(open_file):
        
        if 'D I S P L A C E M E N T   V E C T O R'  in line:

            # skip lines to get to start of data
            for i in range(2):
                
                open_file.readline()
            
            # initialize unsorted displacement field
            u_unsorted = np.zeros([6,n_grids])
            
            # loop through all the grid points
            for j in range(n_grids):
                
                # Split line into multiple parts
                line_split = open_file.readline().split()
                
                # get grid number
                grid = line_split[0]
                
                # find grid number in grids
                index = grids.index(int(grid))
                
                # loop over all 6 degrees of freedom
                for k in range(6):
#                    
                    # get the displacement value
                    u_unsorted[k,index] = line_split[k+2]
                    
            # check grid ordering        
            if grids_order == []:
                
                # use unsorted displacement fields
                u_sorted = u_unsorted
                
            else:
                
                # sort the grids
                u_sorted = sortGridDispl(u_unsorted, grids, grids_order)
            
            # assign displacement field to output
            u[cnt] = u_sorted
            
            # incrememnt index on output vector
            cnt += 1

    open_file.close()

    return u


def importFrequencies(out_dir, out_file, n_modes, n_subcases, debug=True):
    """
    Imports a set of modal frequencies from a F06 file.
    
    Args:
        out_dir: path to the output directory
        out_file: F06 file name
        n_modes: number of modes whose frequencies have to be imported
        debug: flag for debug messages
        
    Returns:
        freq: frequencies (n_modes*1 numpy array)
        
    Remarks:
        1. Modes whose frequencies have to be imported must be sequential.
    """
    
    fname = "pyMSCNastranUtilities.importFrequencies"

    if debug:
        print("Importing frequencies...")
        
    # checking inputs
    
    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")    

    # allocation

    # freq = np.zeros([n_modes,1])

    freq = list(range(n_subcases))

    for i in range(0,n_subcases):

        freq[i] = np.zeros([n_modes])

    # importing frequencies

    open_file = open(os.path.join(out_dir,out_file),"r") 

    while not "R E A L   E I G E N V A L U E S" in open_file.readline():
        pass

    for n in range(0,n_subcases):
        while not 'RADIANS' in open_file.readline():
            pass

        f=np.zeros(n_modes)

        for i in range(0,1):

            line = open_file.readline()

        for i in range(0,n_modes):

            line_split = open_file.readline().split()

            f[i] = line_split[4]

        freq[n] = f

    open_file.close()

    return freq

def sortGridDispl(u_unsorted, grids, grids_order):

    # this function sorts a displacement field according
    # to a user-specified grid order
    # this function can be also used to extract a set of
    # ids from the total imported displacement field 

    # allocation

    u_sorted = np.zeros([len(grids_order),3])

    # loop over the grid ids in the desired order

    for i, grid in enumerate(grids_order):

        # finding position of current grid

        index = grids.index(int(grid))

        # saving displacements

        u_sorted[:,i] = u_unsorted[:,index]

    return u_sorted

def importEigenvectors(out_dir, out_file, n_modes, n_grids, grids,n_subcases,grids_order=[],debug=True):

    # this function imports a set of eigenvectors 
    # from a MSC Nastran f06 output (e.g. SOL103) and in
    # case the user supplies a desired output order does
    # a rearrangement of the values for different grids
    # if no output order is supplied then the values for
    # different grids are returned in the order as from
    # the MSC Nastran output file (ascending order based
    # on grid IDs)

    fname = "pyMSCNastranUtilities.importEigenvectors"

    if debug:
        print("Importing Eigenvectors...")

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.") 

    # allocation
    phi = [[np.zeros([6,n_grids]) for i in range(n_modes)] for j in range(n_subcases)]


    # importing eigenvectors
    open_file = open(os.path.join(out_dir,out_file),"r") 
    
    cnt1 = 0 #  number of eigvectors
    cnt2 = 0 #  number of subcases
    
    for count, line in enumerate(open_file):

                if 'R E A L   E I G E N V E C T O R   N O .'  in line:

                    # skip lines to get to start of data
                    for i in range(2):
                        
                        open_file.readline()
                    
                    # initialize unsorted displacement field
                    u = np.zeros([6,n_grids])
                    
                    # loop through all the grid points
                    for j in range(n_grids):
                        
                        # Split line into multiple parts
                        line_split = open_file.readline().split()
                        
                        # get grid number
                        grid = line_split[0]
                        
                        # find grid number in grids
                        index = grids.index(int(grid))
                        
                        # loop over all 6 degrees of freedom
                        for k in range(6):
                          
                            # get the displacement value
                            u[k,index] = line_split[k+2]


                    # assign displacement field to output
                    phi[cnt2][cnt1] = u
                    
                    # incrememnt index on output vector
                    cnt1+=1
                    
                    if cnt1 == n_modes: # once all mode shapes for each subcase
                                        # are read, move to next subcase
                        
                        cnt1=0
                        cnt2+=1
                        
    open_file.close()
    
    return phi

def importRigidBodyMassData(out_dir, out_file, debug=True):
    """
    Imports the rigid-body mass data from a F06 file. 
    
    Args:
        out_dir: path to the output directory
        out_file: F06 file name
        debug: flag for debug messages
        
    Returns:
        M: FEM model mass (scalar)
        x_G: FEM model center of mass coordinates (3*1 numpy array)
        J_G: FEM model inertia tensor wrt center of mass (3*3 numpy array)
        
    Remarks:
        1. Coordinates are returned in the global coordinate system.
    """

    if debug:
        print("Importing rigid-body mass data...")
        
    M = importMass(out_dir,out_file,debug)
    
    x_G = importCenterOfMass(out_dir,out_file,debug)
    
    J_G = importInertiaTensor(out_dir,out_file,debug)
    
    if debug:
        print("Rigid-body mass data imported")    
    
    return M, x_G, J_G
        

def importMass(out_dir, out_file, debug=True):
    """
    Imports the model mass from a F06 file.
    
    Args:
        out_dir: path to the output directory
        out_file: F06 file name
        debug: flag for debug messages
        
    Returns:
        M: model mass (scalar)
    """
    
    fname = "pyMSCNastranUtilities.importMass"    

    if debug:
        print("Importing total mass...")

    # checking inputs

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")

    open_file = open(os.path.join(out_dir,out_file),"r") 

    while not "MASS AXIS SYSTEM (S)" in open_file.readline():
        pass

    line_split = open_file.readline().split()

    M = float(line_split[1])

    open_file.close()
    
    if debug:
        print("Total mass imported")

    return M


def importCenterOfMass(out_dir, out_file, debug=True):
    """
    Imports the model center of mass coordinates from a F06 file.
    
    Args:
        out_dir: path to the output directory
        out_file: F06 file name
        debug: flag for debug messages
        
    Returns:
        x_G: model center of mass coordinates 
        
    Remarks:
        1. Center of mass coordinates are wrt the GRDPNT reference point.
        2. Coordinates are returned in the global coordinate system.
    """

    fname = "pyMSCNastranUtilities.importCenterOfMass"    

    if debug:
        print("Importing center of mass coordinates...")
        
    # checking inputs

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")        

    # allocation

    x_G = np.zeros([3,1])

    # importing coordinates

    open_file = open(os.path.join(out_dir,out_file),"r") 

    while not "MASS AXIS SYSTEM (S)" in open_file.readline():
        pass

    line_split = open_file.readline().split()

    x_G[1,0] = float(line_split[3])
    x_G[2,0] = float(line_split[4])

    line_split = open_file.readline().split()

    x_G[0,0] = float(line_split[2])

    open_file.close()
    
    if debug:
        print("Center of mass coordinates imported")

    return x_G


def importInertiaTensor(out_dir, out_file, debug=True):
    """
    Imports the model inertia tensor from a F06 file.
    
    Args:
        out_dir: path to the output directory
        out_file: F06 file name
        debug: flag for debug messages
        
    Returns:
        J_G: model inertia tensor 
        
    Remarks:
        1. Inertia tensor is always about the center of mass.
        2. Components are returned in the global coordinate system.
    """
    
    fname = "pyMSCNastranUtilities.importInertiaTensor"    
    
    if debug:
        print("Importing inertia tensor with respect to center of mass...")
        
    # checking inputs

    if not os.path.isfile(os.path.join(out_dir,out_file)):
        raise ValueError("Error in "+fname+": file not found.")        

    # allocation

    J_G = np.zeros([3,3])

    # importing inertia tensor

    open_file = open(os.path.join(out_dir,out_file),'r') 

    while not "I(S)" in open_file.readline():

        pass

    for i in range(0,3):

        line_split = open_file.readline().split()

        for j in range(i,3):

            value = float(line_split[j+1])

            if i == j:

                J_G[i,j] = +value

            else:

                J_G[i,j] = -value
                J_G[j,i] = -value

    open_file.close()
    
    if debug:
        print("Inertia tensor with respetc to center of mass imported")

    return J_G


def writeGRID(out_file, grid_id, grid_coords, cid_inp=0, cid_out=0):
    """
    Writes a GRID entry.
    
    Args:
        out_file: open output file handler 
        grid_id: Grid ID
        grid_coords: grid coordinates
        cid_inp: ID of the coordinate system grid coordinates are given in
        cid_out: ID of the coordinate system grid displacement are output
    """

    # getting grid coordinates

    x = grid_coords[0]
    y = grid_coords[1]
    z = grid_coords[2]

    # writing entry

    out_file.write("$\n")
    out_file.write("GRID,%i,%i,%e,%e,%e,%i\n" %(grid_id,cid_inp,x,y,z,cid_out))

    return


def writeCORD2R(out_file, cid, grid_coords, cid_inp=0):
    """
    Writes a CORD2R entry.
    
    Args:
        out_file: open output file handler 
        cid: coordinate system ID
        grid_coords: grid coordinates defining the coordinate system
        cid_inp: ID of the coordinate system grid coordinates are given in
    """

    # getting coordinates of origin
    
    x1 = grid_coords[0,0] 
    y1 = grid_coords[0,1]
    z1 = grid_coords[0,2]

    # getting coordinates of point along z axis

    x2 = grid_coords[1,0]
    y2 = grid_coords[1,1]
    z2 = grid_coords[1,2]

    # getting coordinates of point in the xz plane

    x3 = grid_coords[2,0]
    y3 = grid_coords[2,1]
    z3 = grid_coords[2,2] 

    # writing entry

    out_file.write("$\n")
    out_file.write("CORD2R,%i,%i,%e,%e,%e,%e,%e,%e\n" %(cid,cid_inp,x1,y1,z1,x2,y2,z2))
    out_file.write(",%e,%e,%e\n" %(x3,y3,z3))

    return


def writeRBE3(out_file, rbe3_id, dep_grid_id, dep_comps, ind_grid_ids, inp_comps):
    """
    Writes a RBE3 entry.
    
    Args:
        out_file: open output file handler 
        rbe3_id: RBE3 ID
        dep_grid_id: dependent grid ID
        dep_comps: dependent degrees of freedom
        ind_grid_ids: independent grid IDs 
        ind_comps: independent degrees of freedom
    """

    # writing entry

    out_file.write("$\n")
    out_file.write("RBE3,%i,,%i,%i,1.,%i" %(rbe3_id,dep_grid_id,dep_comps,inp_comps))

    # initializing field counter

    k = 7
    
    # loop on the independent grids

    for ind_grid_id in ind_grid_ids:

        # incrementing field counter

        k = k+1

        # writing entry continuation

        if not k == 9 and not ind_grid_id == ind_grid_ids[-1]:

            out_file.write(",%i" %ind_grid_id)

        elif k == 9:

            out_file.write(",%i\n" %ind_grid_id)

            k = 1

        elif ind_grid_id == ind_grid_ids[-1]:

            out_file.write(",%i\n" %ind_grid_id)

    return 


def writeRBE2(out_file, rbe2_id, ind_grid_id, dep_grid_ids, dep_comps):
    """
    Writes a RBE2 entry.
    
    Args:
        out_file: open output file handler 
        rbe2_id: RBE2 ID
        ind_grid_id: independent grid ID
        dep_grid_ids: dependent grid IDs
        dep_grid_comps: dependent degrees of freedom
    """

    # writing entry

    out_file.write("$\n")
    out_file.write("RBE2,%i,%i,%i" %(rbe2_id,ind_grid_id,dep_comps))

    # initializing field counter

    k = 4
    
    # loop on the dependent grids

    for dep_grid_id in dep_grid_ids:

        # incrementing field counter

        k = k+1

        # writing entry continuation

        if not k == 9 and not dep_grid_id == dep_grid_ids[-1]:

            out_file.write(",%i" %dep_grid_id)

        elif k == 9:

            out_file.write(",%i\n" %dep_grid_id)

            k = 1

        elif dep_grid_id == dep_grid_ids[-1]:

            out_file.write(",%i\n" %dep_grid_id)

    return


def writeLoad(out_file, load_type, load_id, grid_id, load_comps, cid_inp=0):
    """
    Writes a FORCE or MOMENT entry.
    
    Args:
        out_file: open output file handler 
        load_type: type of entry to be written (FORCE or MOMENT)
        grid_id: grid ID 
        load_comps: load components 
        cid_inp: ID of the coordinate system load is given in
    """
    
    fname = "pyMSCNastranUtilities.writeLoad"

    # checking load type

    if not load_type in ["FORCE", "MOMENT"]:
        raise ValueError("Error in "+fname+": invalid load type.")

    # getting load components

    x = load_comps[0]
    y = load_comps[1]
    z = load_comps[2]

    # writing entry

    out_file.write("$\n")
    out_file.write(load_type+",%i,%i,%i,1.,%e,%e,%e\n" %(load_id,grid_id,cid_inp,x,y,z))

    return

def ComputeMAC(phi_1, phi_2):
    """
    Computes the MAC between two eigenvectors.
    Args:
        phi_1: eigenvector 1 (n_elems*1 numpy array)
        phi_2: eigenvector 2 (n_elems*1 numpy array)
    
    Returns: 
        MAC: MAC between the two eigenvectors
    """
    # check inputs
    if len(phi_1) != len(phi_2):
        raise ValueError("Input eigenvectors have different length.")

    # get number of components
    n_elems = len(phi_1)

    # initialize MAC
    MAC = 0.0

    # get norms of the eigenvectors
    phi_1n = np.linalg.norm(phi_1)
    phi_2n = np.linalg.norm(phi_2)

    # compute MAC
    for i in range(0, n_elems):

        # get elements
        phi_1_i = phi_1[i]
        phi_2_i = phi_2[i]

        # get complex conjugates
        phi_1_ic = np.conj(phi_1_i)

        # add MAC contribution
        MAC = MAC + np.abs(phi_1_ic * phi_2_i)

    MAC = MAC ** 2.0 / (phi_1n * phi_2n)

    return MAC


def ComputeMACs(Phi_1, Phi_2):
    """
   Compute MACs between two eigenvector sets.
    Args:
        Phi_1: eigenvector 1 (n_elems*n_modes numpy array)
        Phi_2: eigenvector 2 (n_elems*n_modes numpy array)
    
    Returns: 
        MACs: MACs between the two eigenvector sets
        
    Remarks:
        1. MAC is computed between the eigenvectors in the two sets that have
           the same order.
    """
    # check inputs
    if len(Phi_1[:, 0]) != len(Phi_2[:, 0]):
        raise ValueError("Eigenvectors have different length.")
    if len(Phi_1[0, :]) != len(Phi_2[0, :]):
        raise ValueError("Eigenvector sets have different size.")

    # get number of modes
    n_modes = len(Phi_1[0, :])

    # allocate MACs
    MACs = np.zeros(n_modes)

    # loop eigenvectors
    for i in range(0, n_modes):

        # compute MAC
        MACs[i] = ComputeMAC(Phi_1[:, i], Phi_2[:, i])

    return MACs


def TrackModes(Phi_1, Phi_2, roots_2):
    """
    Reorders eigenvector sets based on MAC criterion.
    
    Args:
        Phi_1: old eigenvectors (n_elems*n_modes numpy array)
        Phi_2: new eigenvectors (n_elems*n_modes numpy array)
        roots_2: new eigenvalues (n_modes*1 numpy array)
        
    Returns:
        Phi_tracked: reordered eigenvectors from set 2 
        roots_tracked: reordered eigenvalues from set 2
        MAC: vector of MAC between set 2 and set 1 eigenvectors
        best_indices: tracking indices of set 2 with respect to set 1
    """
    # check inputs
    if len(Phi_1[:, 0]) != len(Phi_2[:, 0]):
        raise ValueError("Eigenvectors have different length.")
    if len(Phi_1[0, :]) != len(Phi_2[0, :]):
        raise ValueError("Eigenvector sets have different size.")

    # get number of components
    n_elems = len(Phi_1[:, 0])

    # get number of modes
    n_modes = len(Phi_1[0, :])

    # allocate tracked roots
    roots_tracked = np.zeros(n_modes)

    # allocate tracked eigenvectors
    Phi_tracked = np.zeros([n_elems, n_modes])

    # allocate MAC
    MAC = np.zeros(n_modes)

    # initialize list of matched indices
    best_indices = []

    # loop new eigenvectors
    for i in range(0, n_modes):

        # initialize best MAC
        MAC_ij_best = 0.0

        # loop old eigenvectors
        for j in range(0, n_modes):

            # if not already matched
            if not j in best_indices:

                # compute MAC
                MAC_ij = ComputeMAC(Phi_1[:, j], Phi_2[:, i])

                # update MAC if applicable
                if MAC_ij > MAC_ij_best:

                    MAC_ij_best = MAC_ij
                    MAC_ij_best_index = j

        # add matched index to list to be excluded
        best_indices.append(MAC_ij_best_index)

        # save tracked root
        roots_tracked[MAC_ij_best_index] = roots_2[i]

        # save tracked eigenvector
        Phi_tracked[:, MAC_ij_best_index] = Phi_2[:, i]

        # save MAC
        MAC[MAC_ij_best_index] = MAC_ij_best

    return Phi_tracked, roots_tracked, MAC, best_indices