import numpy as np, sys

# Test assertion function
def assert_passed(init_hsize, hsize, init_V, V, interactive_run=hasattr(sys, "ps1")):
    error_str = "Error in flipCell when trying to flip cell with hSize={:} and volume V={:.2e} (gave a new hSize={:} and V={:.2e})".format(init_hsize, init_V, hsize, V)
    if not interactive_run: raise YadeCheckError(error_str)
    else: print(error_str)

# Test parameters
ucube = Matrix3(np.eye(3)) # unit cube, to be recovered after each flipCell
max_nflips = 5 # maximum number of flips necessary to recover the unit cube
n_random_tests = 1e3 # Number of random flips to be performed
fail_flag = False
O.periodic = True
maxTries=3
multiFlips=0 # count how many cases needed multiple flips

O.engines= [InsertionSortCollider()] # so we don't spam logs with "<WARNING> Cell:95 yade::Matrix3r yade::Cell::flipCell(): No collider found"

# Test, for all directions, flipping in only one direction with both orientations
for ax1 in range(3): 
    for ax2 in range(3): # iterate over all possible flip directions
        if ax1==ax2: continue
        for ori in [-1,1]: # iterate over orientations
            for nflips in range(1, max_nflips+1):
                # Set hSize to a flipped unit cube
                hSize, flip_matrix = np.array(ucube), np.array(ucube)
                hSize[ax1, ax2] = ori*nflips
                flip_matrix[ax1, ax2] = -ori*nflips

                # Test flipCell
                O.cell.hSize = Matrix3(hSize)
                init_V = O.cell.volume
                applied_flip_matrix = flipCell()

                # Check results
                H = O.cell.hSize
                orthogonal = Vector3(H[0]).cross(Vector3(H[1])).cross(Vector3(H[2])).norm() < 1e-10
                if not orthogonal or O.cell.volume!=1: 
                    print("ortho=",Vector3(H[0]).cross(Vector3(H[1])).cross(Vector3(H[2])).norm())
                    fail_flag = True
                    assert_passed(Matrix3(hSize), O.cell.hSize, init_V, O.cell.volume)

# Random tests on several flips along different directions
for test in range(int(n_random_tests)):
    hSize = np.array(ucube)

    # Flip the cell manually
    for flip in range(1, max_nflips+1):
        ax1, ax2, _ = np.random.permutation([0,1,2]) # get two random axes
        ori = np.random.choice([-1, 1]) # get a random orientation
        hSize[:,ax1] += ori * hSize[:,ax2] # flip the cell

    # Set new hSize
    O.cell.hSize = Matrix3(hSize.copy())
    init_V = O.cell.volume
    H = O.cell.hSize
    orthogonal = Vector3(H[0]).cross(Vector3(H[1])).cross(Vector3(H[2])).norm() < 1e-10

    # Test flipCell
    numTries = 0
    while (not orthogonal) and numTries<maxTries:
        applied_flip_matrix = flipCell()
        # Check results
        H = O.cell.hSize
        orthogonal = Vector3(H[0]).cross(Vector3(H[1])).cross(Vector3(H[2])).norm() < 1e-10
        numTries += 1
        
    if numTries>1: multiFlips+=1
    if not orthogonal or O.cell.volume!=1:
        print("ortho=",Vector3(H[0]).cross(Vector3(H[1])).cross(Vector3(H[2])).norm())
        fail_flag = True
        assert_passed(Matrix3(hSize), O.cell.hSize, init_V, O.cell.volume)

if not fail_flag: print("flipCell is working properly (tested on {:d} random flips, {:d} of which needed more than one flip for orthogonality)".format(int(n_random_tests),int(multiFlips)))
 
