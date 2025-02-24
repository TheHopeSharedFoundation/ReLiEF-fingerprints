import numpy as np
import sklearn
from sklearn.decomposition import PCA
from scipy.sparse import coo_matrix
import scoria
import mdtraj as md

traj = md.load('/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.xtc', top='/Users/benjaminsamudio/3w32-benzene_NPT_production_2024sep13utc045116_transformed_skip_0.pdb')
pocket = scoria.Molecule()
platform = scoria.Molecule()
pocket.load_pdb_into("/Users/benjaminsamudio/Desktop/EGFR_inactiveState_pocketOfInterest.pdb")
platform.load_pdb_into("/Users/benjaminsamudio/Desktop/Ne_atom_grid_0p50_angstrom_spacing_2d_testing.pdb")

translation_vector = np.array([0,0,-0.5])
reset_vector = np.array([0,0,25])
for frame in range(1,int(traj.n_frames)+1):
        single_frame = traj[frame]
        single_frame.save_pdb('single_frame.pdb')
        protein = scoria.Molecule()
        protein.load_pdb_into("single_frame.pdb")
        print(f"Operating on frame: {frame}")
        row = []
        column = []
        value = []
        for i in range(0,50):
                platform.translate_molecule(translation_vector)
                platform_clash_protein = platform.select_close_atoms_from_different_molecules(protein,1.8)
                platform_avoid_protein = platform.invert_selection(platform_clash_protein[0])
                new_mol = scoria.Molecule()
                new_mol =  platform.get_molecule_from_selection(platform_avoid_protein)
                close_reference_pocket = new_mol.select_close_atoms_from_different_molecules(pocket,0.6)
                if len(close_reference_pocket[0]) > 0:
                        for j in list(close_reference_pocket[0]):
                                row.append(int(frame))
                                column.append(i*1600+j)
                                value.append(1)
        print(row)
        print(column)
        print(value)
        platform.translate_molecule(reset_vector)


row_array = np.array([3,10,12,15]) #<----- max 16
column_array = np.array([4,10,30,32]) #<------ max 33
data_array = np.array([1,1,1,1])

print(row_array)
print(column_array)
print(data_array)

sparse_matrix = coo_matrix((data_array, (row_array, column_array)), shape=(16, 33))
print(sparse_matrix)
dense_matrix = sparse_matrix.toarray()
pca = PCA(n_components=2)  # Extract the top 2 principal components
transformed_data = pca.fit_transform(dense_matrix)
for i in range(0,len(transformed_data)):
        print(f"{transformed_data[i][0]},{transformed_data[i][1]}")

print(transformed_data)
