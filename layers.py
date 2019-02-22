'''
The layer class deals with the properties of each layer
of a film. The film has a thickness d and a complex refractive
index = n + i*k. If the material is already saved, n and k do not need
to be input as parameters when instantiating the class.
'''
import os
import pandas as pd
from tabulate import tabulate

class layer:
    def __init__(self, name, d, n = 0, k = 0):
        #Setup layer properties
        self.name = name
        self.load_data_from_file(n, k)
        self.d = d

        #Variables will be fit for when set to true
        self.fit_n = False
        self.fit_k = False
        self.fit_d = False
        self.n_guess = None
        self.k_guess = None
        self.d_guess = None

    def load_data_from_file(self, n, k):
        #Set n and k if user inputs them
        if n:
            self.n = n
        if k:
            self.k = k

        #Check for storage file if n and k not input
        if (not (n and k)) and os.path.isfile('./Material_Information.txt'):
            data = pd.read_csv('./Material_Information.txt', sep = '\t', header = 0)
        elif not (self.k and self.n):
            raise Exception('Material_Information.txt not found.')

        #Read data from file if not input
        try:
            data = data[data.Material == self.name]
            if not data.empty:
                if not n:
                    self.n = data.loc[data.index, 'n'].iat[0]
                if not k:
                    self.k = data.loc[data.index, 'k'].iat[0]
            else:
                raise Exception('Material not found in file and no n/k values provided.')
        except:
            return

    def save_material(self):    
        #Load material data
        if os.path.isfile('./Material_Information.txt'):
            data = pd.read_csv('./Material_Information.txt', sep = '\t', header = 0)
        else:
            #Create file, add headers, and load from pandas
            temp_file = open('Material_Information.txt', 'a')
            temp_file.write('Material\tn\tk')
            temp_file.close()
            data = pd.read_csv('./Material_Information.txt', sep = '\t', header = 0)
        
        #Don't add if already in file
        if self.name in list(data['Material']):
            print('Material is already in file.')
            return

        #Create new DataFrame, append to data, and write to file
        string = pd.DataFrame([[str(self.name), str(self.n), str(self.k)]], \
                              columns = ['Material', 'n', 'k'])
        data = data.append(string)
        data.to_csv('Material_Information.txt', sep = '\t', index = False)
        return

    def remove_material(self):    
        #Load material data
        if os.path.isfile('./Material_Information.txt'):
            data = pd.read_csv('./Material_Information.txt', sep = '\t', header = 0)
        else:
            print('File does not exist to remove data from.')
            return
        
        #Don't remove if not in file
        if self.name not in list(data['Material']):
            print('Material is not in file.')
            return

        #Remove data and write to file
        data = data[data.Material != self.name]
        data.to_csv('Material_Information.txt', sep = '\t', index = False)
        return

    def update_material(self):
        #Load material data
        if os.path.isfile('./Material_Information.txt'):
            data = pd.read_csv('./Material_Information.txt', sep = '\t', header = 0)
        else:
            print('Material_Information.txt not found, use save_material to create one.')
            return

        #Find index of material and change n/k values
        index = data[data['Material'] == self.name].index.item()
        data.at[index, 'n'] = self.n
        data.at[index, 'k'] = self.k

        #Write to file
        data.to_csv('Material_Information.txt', sep = '\t', index = False)
        return

    def set_fit_variables(self, properties):
        #Flag property to be fit in solver
        if any('n' in sublist for sublist in properties):
            self.fit_n = True
            self.n_guess = [item[1] for item in properties if item[0] == 'n'][0]
        if any('k' in sublist for sublist in properties):
            self.fit_k = True
            self.k_guess = [item[1] for item in properties if item[0] == 'k'][0]
        if any('d' in sublist for sublist in properties):
            self.fit_d = True
            self.d_guess = [item[1] for item in properties if item[0] == 'd'][0]

    def get_fit_variables(self):
        variables = []

        if self.fit_n:
            variables.append(['n', self.n_guess])
        if self.fit_k:
            variables.append(['k', self.k_guess])
        if self.fit_d:
            variables.append(['d', self.d_guess])

        return variables

    def __str__(self):
        return tabulate([['n', str(self.n), str(self.fit_n), str(self.n_guess)], \
                         ['k', str(self.k), str(self.fit_k), str(self.k_guess)], \
                         ['d', str(self.d), str(self.fit_d), str(self.d_guess)]], \
                         headers = [self.name, 'True', 'Fit?', 'Guess'])
