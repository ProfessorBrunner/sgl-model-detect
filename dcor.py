import math

'''
A Class for calculating the distance correlation beween two vectors of 'same' size.
Constructor	: 	dcor(X,Y)				Returns an object of dcor class. X and Y are input vectors.
Functions  	: 	find_correlation()		Returns the distance correlation.
'''

class dcor:
	def __init__(self, X, Y):
		X_Mat = [] 						# X distance matrix - for all pairwise distances in X 
		Y_Mat = []						# Y distance matrix - for all pairwise distances in Y
		D_X_Mat = []					# X doubly centered distance matrix
		D_Y_Mat = []					# Y doubly centered distance matrix
		X_Row_Sum = []					# i'th Row sum of X distance matrix
		Y_Row_Sum = []					# i'th Row sum of Y distance matrix
		X_Col_Sum = []					# i'th Col sum of X distance matrix
		Y_Col_Sum = []					# i'th Col sum of Y distance matrix
		self.dCor = 0					# class variable - Correlation between X and Y

		if(len(X)!=len(Y)):				# Check for same size of X and Y
			print "Distance Correlation: Input vectors must be of same size -"+str(len(X))+" != "+str(len(Y))
			return -1	
		
		if(len(X)==0):					# If Empty
			print "Distance Correlation: Input vector empty"
			return -1		
		
		for i in range (0,len(X)):
			cur_Row_x = []				# Temporary variables for forming a 2D Matrix
			cur_Row_y = []
			sum_r_x = 0
			sum_r_y = 0			
			for j in range (0,len(X)):
				sum_r_x = sum_r_x + abs( X[i] - X[j])			# Adding sum to find the row mean 
				sum_r_y = sum_r_y + abs( Y[i] - Y[j]) 				
				cur_Row_x.append( abs( X[i] - X[j])  )			# Eucledian norm equivalent to modulus here
				cur_Row_y.append( abs( Y[i] - Y[j])  )
			
			X_Row_Sum.append(sum_r_x) 	
			Y_Row_Sum.append(sum_r_y) 	
			X_Mat.append(cur_Row_x)
			Y_Mat.append(cur_Row_y)

		for i in range (0,len(X)):								# Finding the column mean
			sum_c_x = 0
			sum_c_y = 0			
			for j in range (0, len(X)):
				 sum_c_x += X_Mat[j][i]		
				 sum_c_y += Y_Mat[j][i]
			X_Col_Sum.append(sum_c_x)
			Y_Col_Sum.append(sum_c_y) 			


		X_mean = sum(X_Row_Sum) / (len(X) ** 2.0)					# Overall mean equivalent to sum of all row sums  
		Y_mean = sum(Y_Row_Sum) / (len(X) ** 2.0)					# divided by total number of elements i.e. n^2

		for i in range (0,len(X)):
			cur_row_x = []
			cur_row_y = []
			for j in range (0, len(X)):								# Forming doubly centered distance matrix - cast len(X) to float for div.
				cur_row_x.append( X_Mat[i][j] - X_Row_Sum[i] / float(len(X)) - X_Col_Sum[j] / float(len(X)) + X_mean )
 				cur_row_y.append( Y_Mat[i][j] - Y_Row_Sum[i] / float(len(X)) - Y_Col_Sum[j] / float(len(X)) + Y_mean )
 			D_X_Mat.append(cur_row_x)
			D_Y_Mat.append(cur_row_y)

		Sum_Cov = 0
		Sum_Var_X = 0
		Sum_Var_Y = 0
		for i in range (0,len(X)):
			for j in range (0, len(X)):
				Sum_Cov 	= Sum_Cov 	+ D_X_Mat[i][j] * D_Y_Mat[i][j]
				Sum_Var_X	= Sum_Var_X	+ D_X_Mat[i][j] * D_X_Mat[i][j]
				Sum_Var_Y 	= Sum_Var_Y	+ D_Y_Mat[i][j] * D_Y_Mat[i][j]
		
		dCov 	= Sum_Cov 	/ (len(X) ** 2.0)
 		dVar_X 	= Sum_Var_X / (len(X) ** 2.0)
 		dVar_Y 	= Sum_Var_Y / (len(X) ** 2.0)
 	
		self.dCor = dCov / math.sqrt( dVar_X * dVar_Y)

	def find_correlation(self):
		return self.dCor			 
