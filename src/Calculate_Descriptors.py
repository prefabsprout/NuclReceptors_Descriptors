from COM import COM_protein
from COM_clamp import COM_clamp
from COM_helix import COM_helix
from PairwiseSep import PairwiseSep
from COM_Calpha_angle import COM_Calpha_angle
from Prot_hel_dist import prot_hel_dist
from sse_percentage import sseCalc
from len_of_helix import len_of_hel
from cos_hel import cos_hel
from clamp_dist import ch_clamp_dist
from clamp_angle import ch_clamp_angles
import os
import pandas as pd


if __name__ == '__main__':
    # clamp_resid - residues which form a clamp (you need to get it from literature sources
    # species - fill species_name column
    # prep - fill preparation status of your PDB file
    def calc_all(dir, db_output, clamp_resid, species, prep):
        # Create columns for every descriptor
        cols = ['prot_name', 'species_name', "preparation"]

        com_cols = ['prot_name','COM_x', 'COM_y', 'COM_z']

        angle_cols = ['prot_name']
        for elem in range(0, 12):
            angle_cols.append('Angle_between_COM_and_Calpha_of_hel_' + str(elem + 1))

        comclampdist_cols = ['prot_name','dist_between_COM_clamp_1', 'dist_between_COM_clamp_2', 'dist_between_COM_clamp_3']

        comhel_cols = ['prot_name']
        for elem in range(0, 12):
            comhel_cols.append('COM_helix_x_num_' + str(elem + 1))
            comhel_cols.append('COM_helix_y_num_' + str(elem + 1))
            comhel_cols.append('COM_helix_z_num_' + str(elem + 1))

        pairwise_cols = ['prot_name']
        for i in range(0, 11):
            for j in range(i + 1, 12):
                pairwise_cols.append('Pairwise_sep_between_hel_' + str(i + 1) + '_and_hel_' + str(j + 1))

        protheldist_cols = ['prot_name']
        for elem in range(0, 12):
            protheldist_cols.append('Prot_Helix_' + str(elem + 1) + '_distance')

        cols_sse = ['prot_name', 'Helix', 'Beta bridge', 'Strand', 'Helix-3', 'Helix-5', 'Turn', 'Bend', 'Other']

        cols_len = ['prot_name']
        for elem in range(1, 13):
            cols_len.append(f'len_of_hel_{elem}')

        cols_cos = ['prot_name']
        for i in range(1, 12):
            for j in range(i + 1, 13):
                cols_cos.append(f'cos_between_hel_{i}_and_hel_{j}')

        cols_cl_dist = ['prot_name']
        for elem in range(1, 4):
            cols_cl_dist.append(f'dist_clamp_{elem}')

        cols_cl_angle = ['prot_name']
        for elem in range(1, 3):
            for el in range(elem + 1, 4):
                cols_cl_angle.append(f'clamp_angle_{elem}-{el}')

        # Create empty dataframes for every descriptor
        # Oh boy, that one is definitely not a golden standard of programing
        df = pd.DataFrame(columns=cols)
        df_clamps = pd.DataFrame(columns=comclampdist_cols)
        df_com = pd.DataFrame(columns=com_cols)
        df_alphaangle = pd.DataFrame(columns=angle_cols)
        df_comhels = pd.DataFrame(columns=comhel_cols)
        df_pairseps = pd.DataFrame(columns=pairwise_cols)
        df_prothel = pd.DataFrame(columns=protheldist_cols)
        df_cl_dist = pd.DataFrame(columns=cols_cl_dist)
        df_cl_angles = pd.DataFrame(columns=cols_cl_angle)
        df_sse = pd.DataFrame(columns=cols_sse)
        df_len = pd.DataFrame(columns=cols_len)
        df_cos = pd.DataFrame(columns=cols_cos)

        for filename in os.listdir(dir):
            data = [filename, species, prep]
            df = df.append(pd.Series(data, index=cols[0:len(data)]), ignore_index=True)

            # Fill every sub-dataframe with data
            try:
                COM_prot = COM_protein(os.path.join(dir, filename))
            except:
                KeyError
            data_com = [filename]
            for coord in COM_prot:
                data_com.append(coord)
            df_com = df_com.append(pd.Series(data_com, index=com_cols[0:len(data_com)]), ignore_index=True)

            try:
                alpha_angle = COM_Calpha_angle(os.path.join(dir, filename))
            except:
                KeyError
            data_alphaagnle = [filename]
            for elem in alpha_angle:
                data_alphaagnle.append(elem)
            df_alphaangle = df_alphaangle.append(pd.Series(data_alphaagnle, index=angle_cols[0:len(data_alphaagnle)]), ignore_index=True)

            try:
                clamps = COM_clamp(os.path.join(dir, filename), clamp_resid)
            except:
                KeyError
            data_clamps = [filename]
            if clamps is not None:
                for dist in clamps:
                    data_clamps.append(dist)
            df_clamps = df_clamps.append(pd.Series(data_clamps, index=comclampdist_cols[0:len(data_clamps)]), ignore_index=True)

            try:
                com_hels = COM_helix(os.path.join(dir, filename))
            except:
                KeyError
            data_comhels = [filename]
            for elem in com_hels:
                for coord in elem[0]:
                    data_comhels.append(coord)
            df_comhels = df_comhels.append(pd.Series(data_comhels, index=comhel_cols[0:len(data_comhels)]), ignore_index=True)

            try:
                pairseps = PairwiseSep(os.path.join(dir, filename))
            except:
                KeyError
            data_pairseps = [filename]
            for elem in pairseps:
                for dist in elem:
                    data_pairseps.append(dist)
            df_pairseps = df_pairseps.append(pd.Series(data_pairseps, index=pairwise_cols[0:len(data_pairseps)]), ignore_index=True)

            try:
                prothel = prot_hel_dist(os.path.join(dir, filename))
            except:
                KeyError
            data_prothel = [filename]
            for elem in prothel:
                data_prothel.append(elem)
            df_prothel = df_prothel.append(pd.Series(data_prothel, index=protheldist_cols[0:len(data_prothel)]), ignore_index=True)

            sse = None
            try:
                sse = sseCalc(os.path.join(dir, filename))
            except KeyError:
                pass
            data_sse = [filename]
            if sse is not None:
                for struct in sse:
                    data_sse.append(sse[struct])
            df_sse = df_sse.append(pd.Series(data_sse, index=cols_sse[0:len(data_sse)]), ignore_index=True)

            lens_hels = None
            try:
                lens_hels = len_of_hel(os.path.join(dir, filename))
            except KeyError:
                pass
            data_lens = [filename]
            if lens_hels is not None:
                for lens in lens_hels:
                    data_lens.append(lens_hels[lens])
            df_len = df_len.append(pd.Series(data_lens, index=cols_len[0:len(data_lens)]), ignore_index=True)

            cos = None
            try:
                cos = cos_hel(os.path.join(dir, filename))
            except KeyError:
                pass
            data_cos = [filename]
            if cos is not None:
                for elem in cos:
                    data_cos.append(cos[elem])
            df_cos = df_cos.append(pd.Series(data_cos, index=cols_cos[0:len(data_cos)]), ignore_index=True)

            clamp_dist = None
            try:
                clamp_dist = ch_clamp_dist(os.path.join(dir, filename), clamp_resid)
            except KeyError:
                pass
            cl_dist = [filename]
            if clamp_dist is not None:
                for elem in clamp_dist:
                    cl_dist.append(clamp_dist[elem])
            df_cl_dist = df_cl_dist.append(pd.Series(cl_dist, index=cols_cl_dist[0:len(cl_dist)]),
                                           ignore_index=True)

            clamp_angle = None
            try:
                clamp_angle = ch_clamp_angles(os.path.join(dir, filename), clamp_resid)
            except KeyError:
                pass
            cl_angle = [filename]
            if clamp_angle is not None:
                for elem in clamp_angle:
                    cl_angle.append(clamp_angle[elem])
            df_cl_angles = df_cl_angles.append(pd.Series(cl_angle, index=cols_cl_angle[0:len(cl_angle)]),
                                               ignore_index=True)


            # And concatenate it all
            df_concat = df.merge(df_prothel, on='prot_name').merge(df_pairseps,on='prot_name')\
                .merge(df_comhels,on='prot_name').merge(df_clamps,on='prot_name')\
                .merge(df_alphaangle,on='prot_name').merge(df_com,on='prot_name')\
                .merge(df_sse, on='prot_name').merge(df_len, on='prot_name')\
                .merge(df_cos, on='prot_name').merge(df_cl_dist, on='prot_name')\
                .merge(df_cl_angles, on='prot_name')

        df_concat.to_csv(db_output)


    calc_all("raw_data/VDR_PDB/Danio rerio", "calc_results_danio.csv", [274, 292, 446],
             "zebrafish",0)
    calc_all("raw_data/VDR_PDB/Homo sapiens", "calc_results_homosap.csv", [246, 264, 420],
             "human",0)
    calc_all("raw_data/VDR_PDB/Rattus norvegicus", "calc_results_rattus.csv", [242, 260, 416],
             "rat",0)

    calc_all("raw_data/VDR_PDB_prep/Danio rerio", "calc_results_danio_prep.csv",[274, 292, 446],
             "zebrafish",1)
    calc_all("raw_data/VDR_PDB_prep/Homo sapiens", "calc_results_homosap_prep.csv",[246, 264, 420],
             "human",1)
    calc_all("raw_data/VDR_PDB_prep/Rattus norvegicus", "calc_results_rattus_prep.csv",[242, 260, 416],
             "rat",1)
