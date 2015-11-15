#Caluclations for input parameters for cavety flow solutions

#INPUTS
n = 129			#Grid size
Re = 100		#Reynolds Number
d = 1			#Wall height
U = 1			#Freestream velocity


nu = d * U / Re 	#dynamic viscosity
print 'Dynamic Viscosity = ', nu
print nu
dx = d / n
print 'dx = ', dx
dt = 1 * dx / U		#from CFL number
print 'dt = ', dt

