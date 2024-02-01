import numpy as np

def n_glass(x):
    n=(1+0.614555251/(1-0.0145987884/x**2)+0.656775017/(1-0.00287769588/x**2)+1.02699346/(1-107.653051/x**2))**.5
    return n

def n_air(x):
    n=1+0.05792105/(238.0185-x**-2)+0.00167917/(57.362-x**-2)
    return n

def n_water(x): 
    n  = 1.3317
    return n

def parabola(x, a, b, c):
    return a * x**2 + b * x + c

class Index:
    '''
    class to obtain a certain index(indices) of a given vector 
    '''

    def __init__(self,vector):
        self.vector = vector
        return None
    
    def get_index(self, value):
        '''
        method to get the nearest index to a given value
        '''
        absolute_diff = np.abs(self.vector - value)
        index = np.argmin(absolute_diff)

        return index

    def get_double_index(self, value):
        ''' 
        method to get the 2 nearest indices to a given value
        (intended to be used when the function is bivaluated)
        '''
        absolute_diff = np.abs(self.vector-value)
        index_1 = np.argmin(absolute_diff)
        absolute_diff[index_1-5:index_1+5] = 100
        index_2 = np.argmin(absolute_diff)

        if index_1 > index_2:
            i_1 = index_2
            i_2 = index_1
        else:
            i_1 = index_1
            i_2 = index_2

        return i_1, i_2
    
    def get_triple_index(self, value):
        ''' 
        method to get the 2 nearest indices to a given value
        (intended to be used when the function is bivaluated)
        '''
        absolute_diff = np.abs(self.vector-value)
        index_1 = np.argmin(absolute_diff)
        absolute_diff[index_1-7:index_1+7] = 100
        index_2 = np.argmin(absolute_diff)
        absolute_diff[index_2-7:index_2+7] = 100
        index_3 = np.argmin(absolute_diff)

        sorted_indices = np.sort([index_1, index_2, index_3])

        # Assign the sorted indices
        i_1, i_2, i_3 = sorted_indices

        return i_1, i_2, i_3
    
    def get_delimited_index(self, value, index):

        total_vector = self.vector
        delimited_vector = total_vector[index-15 : index + 15]

        absolute_diff = np.abs(delimited_vector-value)

        index_1 = np.argmin(absolute_diff)

        final_index = index - 15 + index_1

        return final_index



  
    