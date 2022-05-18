# -*- encoding=utf-8 -*-
# 2022 Sacha Duverger <sacha.duverger~a~inrae.fr>
"Demonstrate several usage of the flipCell function"
import numpy as np

# A function to format and print the results of all steps in this example, not relevent to flipCell itself
step, planes = 1, [(0,1), (0,2), (1,2)]
def print_hSize(step_print=True, hSize_types="int"):
	global step
	if step_print: print("\n__________________________________________________________________________________________\nStep {:d}".format(step))

	print("\tCell's base vectors (O.cell.hSize) are:")
	if hSize_types=="int":
		for axis in range(3): print("\t\tAxis {:d}: ({:.0f}, {:.0f}, {:.0f})".format(axis, *O.cell.hSize.col(axis)))
	elif hSize_types=="float":
		for axis in range(3): print("\t\tAxis {:d}: ({:.2e}, {:.2e}, {:.2e})".format(axis, *O.cell.hSize.col(axis)))

	print("\n\tAngles between the cell's base vectors and the global axes on the (i,j) plane")
	for ax1, ax2 in planes: 
		phi_1a = np.arccos(O.cell.hSize[ax1,ax1] / np.sqrt((O.cell.hSize[ax1,ax1]**2 + O.cell.hSize[ax2,ax1]**2)))*180/np.pi
		phi_1b = np.arccos(O.cell.hSize[ax2,ax2] / np.sqrt((O.cell.hSize[ax2,ax2]**2 + O.cell.hSize[ax1,ax2]**2)))*180/np.pi
		print("\t\t(i={:d}, j={:d})".format(ax1, ax2))
		for ax, phi_1 in zip([ax1, ax2], [phi_1a, phi_1b]): print("\t\t\tfor cell's base vector {:d} = {:.2f} degrees".format(ax, phi_1))
	if step_print: step += 1

# Cell flipping is performed on periodic cells, so the simulation has to be periodic
O.periodic = True

# Step 1
# Initially, we setup the cell as an unit cube
O.cell.hSize = Matrix3(1,0,0, 0,1,0, 0,0,1)
print_hSize()

# Step 2
# We can use flipCell to shift the top face (whose outer unit normal is (0,1,0)) of the cell 1 grid point to the left along the axis 0
flipCell(Matrix3(0,-1,0, 0,0,0, 0,0,0))
print_hSize()

# Step 3
# If we call flipCell without input argument, the best flip is automatically sought among the neighbour grid points
# The best flip is the one that minimizes the angles between the cell's base vectors and the global axes
flipCell()
print_hSize() # The unit cube is recovered

# Step 4
# When passing a flip matrix, one can shift the faces of the cell by several grid points (in our case the result is an extremly sheared cell)
flipCell(Matrix3(0,-5,0, 0,0,0, 0,0,0))
print_hSize()

# Step 5
# Several flips can be performed in different directions
flipCell(Matrix3(0,0,1, 0,0,0, 0,-2,0))
print_hSize()

# Step 6
# flipCell returns the flip matrix it used on the cell. Since it only flip the cell on grid node at the time, one might need to call it sevral times
# Also, hSize changes between each call of flipCell, the sum of all flip matrices is not the "overall" flip matrix.
while flipCell()!=Matrix3(0,0,0, 0,0,0, 0,0,0): pass # Note that flipCell performs at least 11 tests when called
print_hSize() # The unit cube is recovered

# Step 7
# flipCell works for any hSize
O.cell.hSize = Matrix3(0.001767398585873738732,-3.402651098340004036e-22,-2.455726386595863492e-23, 2.840638172828420518e-23,0.004990738513508447341,-0.0009967331209356602214, -1.550207193852494163e-22,-0.001442690301487606743,0.00144229021430840711)
print_hSize(hSize_types="float")
flipCell()
print("\nStep 7 - After flipping\n")
print_hSize(False, "float") # The best flip was found

