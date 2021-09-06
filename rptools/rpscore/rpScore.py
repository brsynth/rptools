from os import(
    path as os_path,
    remove as os_remove
)
from sys import exit as sys_exit
from pandas import (
    DataFrame as pd_DataFrame,
    read_csv as pd_read_csv
)
import numpy as np
from rdkit import (
    Chem,
    DataStructs,
    RDLogger
)
from rdkit.Chem import AllChem
from h5py import File as h5py_File
from tqdm import tqdm
from sklearn.utils import shuffle
from statistics import (
    mean,
    stdev
)
from xgboost import DMatrix
from pickle import load as pickle_load

##############################################################################################

def feature_template_df(no_of_rxns_thres):
    '''Creates a template dataframe
       containing all the feature columns'''
    feature_data_columns = []
    for n in range(no_of_rxns_thres):
        smarts = "Rxn"+str(n+1)+"_SMARTS"
        rxn_delta_g = "Rxn"+str(n+1)+"_DeltaG"
        rxn_rule_score = "Rxn"+str(n+1)+"_Rule_Score"
        feature_data_columns.extend([smarts, rxn_delta_g, rxn_rule_score])
    feature_data_columns.extend(["Pathway_Delta_G", "Pathway_Flux", "Pathway_Score", "Round1"])
    #print(feature_data_columns)
    feature_data = pd_DataFrame(  columns = feature_data_columns , index = None)
    return feature_data

def loop(i , temp, data):
    '''Returns the indices of all the reactions
       for a pathway in the dataset'''
    temp_list = []
    break_index = None
    flag = True
    j = 1
    for index in range(i, len(data)):
        if temp == data.loc[index,'Pathway Name'] and data.loc[index,'Reaction'] == "RP"+str(j):
            j =j+1
            temp_list.append(index)
            if index+1 == len(data):
                flag = False
        else:
            break_index = index
            break
    return temp_list , break_index , flag


def pathways_index_list(data):
    '''Returns the indices of all the reactions
       for each of the pathways in dataset'''
    pathways_index_list = []
    i = 0
    flag = True
    while flag:
        temp = data.loc[i,'Pathway Name']
        temp_list, i , flag = loop( i , temp, data)
        pathways_index_list.append(temp_list)

    return pathways_index_list


def transform_into_pathway_features(data, scores, flag, no_of_rxns_thres):
    '''Generates the dataframe containing
       all the features.
       data and scores ate the 2 inputs files
       The reactions are represented in SMILES'''

    df = feature_template_df(no_of_rxns_thres)
    pathways_list = pathways_index_list(data)
    #print(pathways_list)
    #print(len(pathways_list))
    drop_list = []
    print("Transforming into pathway features...")
    for count_p, rxn_list in tqdm(enumerate(pathways_list)) :
        if len(rxn_list) > 10:
            drop_list.append(count_p)
            continue
        # At the level of each reaction reading data file
        for n , index in enumerate(rxn_list):
            #print(index)
            smarts = "Rxn"+str(n+1)+"_SMARTS"
            rxn_delta_g = "Rxn"+str(n+1)+"_DeltaG"
            rxn_rule_score = "Rxn"+str(n+1)+"_Rule_Score"
            df.loc[count_p, smarts] = data.loc[index,'Reaction Rule']
            df.loc[count_p, rxn_delta_g] = data.loc[index,'Normalised dfG_prime_m']
            df.loc[count_p, rxn_rule_score] = data.loc[index,'Rule Score']
        # At the level of pathway reading scores file
        df.loc[count_p, "Pathway_Delta_G"] = scores.loc[count_p, 'dfG_prime_m']
        df.loc[count_p, "Pathway_Flux"] = float( scores.loc[count_p,'FBA Flux'].split(';')[1] )
        df.loc[count_p, "Pathway_Score"] = scores.loc[count_p,'Global Score']
        df.loc[count_p, "Lit"] = scores.loc[count_p,'Lit']
        df.loc[count_p, "Round1"] = scores.loc[count_p,'Round1']
    df = df.drop(drop_list)
    df = df.fillna(0)
    if flag:
        df = df[~(df.Round1 < 0)]
        df["Round1"][df['Round1'] > 0] = 1
        #df.to_csv("raw_df.csv")
        df["Round1_OR"] = df["Round1"]
        df = shuffle(df, random_state = 42).reset_index(drop =True)
        for row in range(len(df)):
            if  df.loc[row , "Lit"] == 1 :
                df.loc[row , "Round1_OR"] = 1
        #df.to_csv("raw_df_csv", index = None)
        #print(df)
    else :
        df["Round1_OR"] = df["Round1"]
    return df

def features_encoding (df, flag, data_predict_file):
    '''Creates a HDF5 file containing
       all the features
       Rnx features are encoded in fingerprints'''
    no_of_rxns = 10
    fp_len = 4096
    rxn_len = fp_len + 2
    pathway_len = 3
    y_len = 1

    if flag == "train":
        sys_exit('Encoding feature for training data not available file data_train.h5 must be present in models folder')
    elif flag == "predict":
        path = data_predict_file
        print("Encodining features for the Test set......")
    if os_path.exists(path):
        os_remove(path)
    f=h5py_File(path, "w")
    dset = f.create_dataset('data', (  0, (rxn_len*no_of_rxns + pathway_len + y_len)),dtype='i2',maxshape=(None,(rxn_len*no_of_rxns + pathway_len + y_len)), compression='gzip')

    for row in tqdm(range(len(df))):
        pathway_rxns = np.array([]).reshape(0, rxn_len * no_of_rxns)
        rxns_list = []
        for rxn_no_ in range(no_of_rxns):
            
            rxn_smiles_index = rxn_no_ * 3
            rxn_dg_index = (rxn_no_ + 1)* 3 -2
            rxn_rule_score_index = (rxn_no_ + 1)* 3 - 1
        
            if  str(df.iloc[row , rxn_smiles_index]) != '0':
                #print(df.iloc[row , rxn_smiles_index])
                rxn_smiles = df.iloc[row , rxn_smiles_index]
                rxn_smiles_list = rxn_smiles.split(">>")
                #print(len(rxn_smiles_list))

                if len(rxn_smiles_list) == 2:

                    sub_smiles = rxn_smiles_list[0]
                    sub_m= Chem.MolFromSmiles(sub_smiles)
                    #print(m)
                    sub_fp = AllChem.GetMorganFingerprintAsBitVect(sub_m, 2, nBits = 2048)
                    sub_arr = np.array([])
                    DataStructs.ConvertToNumpyArray(sub_fp, sub_arr)
                    sub_fp= sub_arr.reshape(1,-1)

                    pro_smiles = rxn_smiles_list[1]
                    pro_m= Chem.MolFromSmiles(pro_smiles)
                    #print(m)
                    pro_fp = AllChem.GetMorganFingerprintAsBitVect(pro_m, 2, nBits = 2048)
                    pro_arr = np.zeros((1,))
                    DataStructs.ConvertToNumpyArray(pro_fp, pro_arr)
                    pro_fp= pro_arr.reshape(1,-1)
                    rxn_fp = np.concatenate([sub_fp , pro_fp]).reshape(1, -1)

                elif len(rxn_smiles_list) < 2:
                     
                    pro_smiles = rxn_smiles_list[0]
                    #print(pro_smiles)
                    pro_m= Chem.MolFromSmiles(pro_smiles)
                    #print(pro_m)
                    pro_fp = AllChem.GetMorganFingerprintAsBitVect(pro_m, 2, nBits = fp_len) # JLF: not good !!
                    pro_arr = np.zeros((1,))
                    DataStructs.ConvertToNumpyArray(pro_fp, pro_arr)
                    rxn_fp= pro_arr.reshape(1,-1)
                else:
                    print("There is a problem with the number of components in the reaction")

            else:
                rxn_fp = np.zeros(fp_len).reshape(1,-1)

            rxn_dg = df.iloc[row , rxn_dg_index].reshape(1,-1)
            rxn_rule_score = df.iloc[row , rxn_rule_score_index].reshape(1,-1)
            rxns_list.extend([rxn_fp, rxn_dg, rxn_rule_score])
            #print(rxn_rule_score)

        pathway_rxns = np.concatenate(rxns_list , axis = 1).reshape(1,-1)
        pathway_dg = df.loc[row, "Pathway_Delta_G"].reshape(1,-1)
        pathway_flux = df.loc[row, "Pathway_Flux"].reshape(1,-1)
        pathway_score = df.loc[row, "Pathway_Score"].reshape(1,-1)
        pathway_y = df.loc[row, "Round1_OR"].reshape(1,-1)
        feature = np.concatenate((pathway_rxns, pathway_dg, pathway_flux, pathway_score, pathway_y), axis =1)
        dset.resize(dset.shape[0]+feature.shape[0], axis=0)
        dset[-feature.shape[0]:]= feature
        #print(pathway_flux)

    return dset


def transform_to_matrix(dset, model_file):
    ''''Transforms the prediction dataset into
        an appropriate matrix which is suitable for
        XGBoost'''
    X_test = dset[:,:-1]
    Y_test = dset[:, -1]
    
    if not os_path.exists(model_file):
        sys_exit(f'{model_file} not found')
    else:
        trained_model = pickle_load(open(model_file,'rb'))

    dset_matrix = DMatrix(X_test, label = Y_test)
    trained_model_score =  trained_model.predict(dset_matrix)
    trained_model_score_1 = trained_model_score[:, 1].reshape(-1,1)
    X_test = np.concatenate((X_test, trained_model_score_1), axis = 1)
    dset_matrix = DMatrix(X_test)

    return  dset_matrix
    

###############################################################
def score_prediction(features_dset_train, features_dset_pred, models_path):

    stdev_ = []
    mean_ = []
    pb1_mean = []
    pb1_stdev = []
    print(features_dset_pred)
    n_predictions = np.array([[]])
    print("Predicting n times...")
    for model_number,  n in tqdm( enumerate([ 0, 10, 20, 30, 40, 50, 60, 70,80, 90])): ##########################
        modelfile = os_path.join(
            models_path,
            f'model{model_number}.pickle'
        )
        if not os_path.exists(modelfile):
            sys_exit('modelfile not found')
        model = pickle_load( open(modelfile, 'rb')) ############
        df_test_matrix = transform_to_matrix(
            features_dset_pred,
            os_path.join(
                models_path,
                'model.pickle'
            )
        )
        prediction = model.predict(df_test_matrix)
        pb_1 = prediction[:, 1].reshape(-1, 1)
        prediction = np.asarray([np.argmax(line) for line in prediction]).reshape(-1, 1)
        if n_predictions.shape[1] == 0 :
            n_predictions = prediction
            n_pb_1 = pb_1
        else:
            n_predictions = np.concatenate((n_predictions , prediction), axis = 1)
            n_pb_1 = np.concatenate((n_pb_1, pb_1), axis = 1)

    for row in range(len(n_predictions)):
        line = n_predictions[row, :].tolist()
        line_pb1 = n_pb_1[row, :].tolist()

        mean_.append(mean(line))
        stdev_.append(stdev(line))
        pb1_mean.append(mean(line_pb1))
        pb1_stdev.append(stdev(line_pb1))

    mean_stdev = pd_DataFrame( { 'mean' : mean_ , 'stdev' : stdev_ , 'Prob1_mean': pb1_mean , 'Prob1_stdev' : pb1_stdev})
    mean_stdev.to_csv("mean_stdev_.csv")


##############################################################

# Loading training data saved in models folder
def load_training_data(filename: str):
    if not os_path.exists(filename):
        sys_exit(f'{filename} not found')
    f = h5py_File(filename, "r")
    features_dset_train = f["data"]
    f.close()
    return features_dset_train

# Loading query
def loading_query(
    test_data_file: str,
    test_score_file: str
):
    data_test = pd_read_csv(test_data_file)
    scores_test = pd_read_csv(test_score_file)
    print("Number of pathways : ",len(scores_test))
    print("Total number of reactions  : ", len(data_test))
    return data_test, scores_test

# Encoding and prediction
def encode_and_predict(
    data_test: str,
    scores_test: str,
    data_predict_file: str,
    models_path: str,
    features_dset_train,
    no_of_rxns_thres: int
):
    df_test = transform_into_pathway_features(data_test, scores_test, False, no_of_rxns_thres)
    features_dset_pred  = features_encoding(df_test, "predict", data_predict_file)
    score_prediction(features_dset_train, features_dset_pred, models_path)
    print("Mean - Stdev stats is saved")

def predict_score(
      test_data_file: str,
      test_score_file: str,
      data_predict_file: str,
      models_path: str,
      features_dset_train,
      no_of_rxns_thres: int
):
    data_test, scores_test = loading_query(
      test_data_file,
      test_score_file
    )
    encode_and_predict(
      data_test,
      scores_test,
      data_predict_file,
      models_path,
      features_dset_train,
      no_of_rxns_thres
    )
################################################################


