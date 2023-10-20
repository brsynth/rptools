from os import(
    path as os_path,
    remove as remove
)
from sys import exit as sys_exit
from tempfile import NamedTemporaryFile
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
from logging import (
    Logger,
    getLogger
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
from rptools.rplibs import rpPathway


__CURRENT_PATH = os_path.dirname(
    os_path.abspath(__file__)
)
__MODELS_PATH = os_path.join(
    __CURRENT_PATH,
    'models'
)
__DATA_TRAIN_FILE = os_path.join(
    __MODELS_PATH,
    'data_train.h5'
)


class ThermoError(Exception):
    pass


class FBAError(Exception):
    pass


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


def transform_into_pathway_features(
    data,
    scores,
    flag,
    no_of_rxns_thres
):
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
        df.loc[count_p, "Pathway_Delta_G"] = scores.loc [count_p, 'dfG_prime_m']
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

def features_encoding (df, flag):
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
    # elif flag == "predict":

    print("Encodining features for the Test set......")
    # temp_f_name = str(uuid4())
    f = h5py_File(
        NamedTemporaryFile(delete=True),
        'w'
    )
    number = rxn_len * no_of_rxns + pathway_len + y_len
    dset = f.create_dataset(
        'data',
        (0, number),
        dtype='i2',
        maxshape=(None, number),
        compression='gzip'
    )

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
def score_prediction(
    features_dset_train,
    features_dset_pred,
    models_path,
    out_filename
):

    stdev_ = []
    mean_ = []
    pb1_mean = []
    pb1_stdev = []
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

    mean_stdev = pd_DataFrame( { 'mean': mean_ , 'stdev': stdev_ , 'Prob1_mean': pb1_mean , 'Prob1_stdev': pb1_stdev})
    mean_stdev.to_csv(out_filename)


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
    print("Number of pathways:", len(scores_test))
    print("Total number of reactions:", len(data_test))
    return data_test, scores_test

# Encoding and prediction
def encode_and_predict(
    data_test: str,
    scores_test: str,
    models_path: str,
    features_dset_train,
    no_of_rxns_thres: int,
    out_filename: str
):
    df_test = transform_into_pathway_features(
        data_test,
        scores_test,
        False,
        no_of_rxns_thres
    )
    features_dset_pred  = features_encoding(df_test, "predict")
    score_prediction(
        features_dset_train,
        features_dset_pred,
        models_path,
        out_filename
    )
    # print("Mean - Stdev stats is saved")

def _predict_score(
      test_data_file: str,
      test_score_file: str,
      models_path: str,
      features_dset_train,
      no_of_rxns_thres: int
) -> float:
    # ttdf = open(test_data_file, 'r')
    # print('test_data_file')
    # print(ttdf.read())
    # ttsf = open(test_score_file, 'r')
    # print('test_score_file')
    # print(ttsf.read())

    data_test, scores_test = loading_query(
      test_data_file,
      test_score_file
    )

    with NamedTemporaryFile(delete=False) as out_f:
        encode_and_predict(
            data_test,
            scores_test,
            models_path,
            features_dset_train,
            no_of_rxns_thres,
            out_f.name
        )
        out_f.close()
        score_df = pd_read_csv(out_f.name)
        remove(out_f.name)

    return list(score_df.to_dict()['Prob1_mean'].values())

def predict_score(
    # pathways: List[rpPathway],
    pathway: rpPathway,
    # data_train_file: str,
    # models_path: str,
    no_of_rxns_thres: int,
    logger: Logger = getLogger(__name__)
) -> float:
# ) -> List[float]:

    with NamedTemporaryFile(mode='w', delete=False) as rf:
        with NamedTemporaryFile(mode='w', delete=False) as pf:
            format_files(
                pathway=pathway,
                reactions_fp=rf,
                pathways_fp=pf,
                logger=logger
            )
            rf.close()
            pf.close()
            features_dset_train = load_training_data(
                __DATA_TRAIN_FILE
            )
            score = _predict_score(
                rf.name,
                pf.name,
                __MODELS_PATH,
                features_dset_train,
                no_of_rxns_thres
            )[0]
            remove(rf.name)
            remove(pf.name)
            return score

def format_files(
#   pathways: List[rpPathway],
    pathway: rpPathway,
    reactions_fp,
    pathways_fp,
    logger: Logger = getLogger(__name__)
) -> None:

    # PATHWAYS FILE
    columns = [
        'Unnamed: 0',
        'target_name',
        'target_structure',
        'Pathway Name',
        'pathway_dummy_name',
        'chassis',
        'Reaction',
        'Global Score',
        'Rule ID',
        'EC_number',
        'Reaction Rule',
        'Rule Score',
        'dfG_prime_o',
        'dfG_prime_m',
        'dfG_uncert',
        'Normalised dfG_prime_o',
        'Normalised dfG_prime_m',
        'Normalised dfG_uncert',
        'FBA',
        'FBA Flux',
        'FBA Normalised Flux',
        'UniProt',
        'Selenzyme Score',
        'Measured',
        'Sub',
        'Lit',
        'Round1',
        'Lit_score',
        'Lit_best'
    ]
    # Header
    pathways_fp.write(','.join(columns)+'\n')
    # Body
    i = 0
#   for pathway in pathways:
    data = {}
    data['Unnamed: 0'] = i
    data['Pathway Name'] = pathway.get_id()
    if pathway.get_thermo_dGm_prime():
        data['dfG_prime_m'] = pathway.get_thermo_dGm_prime()['value']
    else:
        raise ThermoError('No thermodynamics data available')
    if pathway.get_fba():
        data['FBA'] = ';'.join(pathway.get_fba().keys())
        data['FBA Flux'] = ';'.join([str(v['value']) for v in pathway.get_fba().values()])
    else:
        raise FBAError('No FBA data available')
    pathways_fp.write(
        ','.join(
        [str(data.get(c, '')) for c in columns]
        ) + '\n'
    )
    i += 1

    # REACTIONS FILE
    columns = [
        '',
        'target_name',
        'target_structure',
        'Pathway Name',
        'pathway_dummy_name',
        'Reaction',
        'Global Score',
        'Rule ID',
        'EC_number',
        'Reaction Rule',
        'Rule Score',
        'dfG_prime_o',
        'dfG_prime_m',
        'dfG_uncert',
        'Normalised dfG_prime_o',
        'Normalised dfG_prime_m',
        'Normalised dfG_uncert',
        'FBA',
        'FBA Flux',
        'FBA Normalised Flux',
        'UniProt',
        'Selenzyme Score',
        'Measured',
        'Sub'
    ]
    # Header
    reactions_fp.write(','.join(columns)+'\n')
    # Body
    i = 0
    # for pathway in pathways:
    # List of reactions with the one producing the target in fist place
    rxn_lst = (
    [pathway.get_target_rxn_id()] +
    [
        rxn_id
        for rxn_id in pathway.get_reactions_ids()
        if rxn_id != pathway.get_target_rxn_id()
    ]
    )
    # smiles = [
    #   '[H]OC(=NC([H])([H])C([H])([H])SC(=O)C([H])=C([H])c1c([H])c([H])c(O[H])c(O[H])c1[H])C([H])([H])C([H])([H])N=C(O[H])C([H])(O[H])C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C1([H])OC([H])(n2c([H])nc3c(N([H])[H])nc([H])nc32)C([H])(O[H])C1([H])OP(=O)(O[H])O[H].[H]OC(=O)C([H])(O[H])C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H]>>CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(O)=NCCC(O)=NCCS.[H]OC(=O)C([H])(OC(=O)C([H])=C([H])c1c([H])c([H])c(O[H])c(O[H])c1[H])C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H]',
    #   '[H]OC(=O)C(=O)C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H].[H+].[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]>>[H]OC(=O)C([H])(O[H])C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H].N=C(O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
    #   '[H]Oc1c([H])c([H])c(C([H])=O)c([H])c1O[H].[H]OC(=NC([H])([H])C([H])([H])SC(=O)C([H])([H])[H])C([H])([H])C([H])([H])N=C(O[H])C([H])(O[H])C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C1([H])OC([H])(n2c([H])nc3c(N([H])[H])nc([H])nc32)C([H])(O[H])C1([H])OP(=O)(O[H])O[H]>>[H]OC(=NC([H])([H])C([H])([H])SC(=O)C([H])=C([H])c1c([H])c([H])c(O[H])c(O[H])c1[H])C([H])([H])C([H])([H])N=C(O[H])C([H])(O[H])C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C1([H])OC([H])(n2c([H])nc3c(N([H])[H])nc([H])nc32)C([H])(O[H])C1([H])OP(=O)(O[H])O[H].[H]O[H]',
    #   '[H]OC(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H].O=O>>[H]OC(=O)C(=O)C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H].N.N',
    #   '[H]Oc1c([H])c([H])c(C([H])=O)c([H])c1[H].O=O>>[H]Oc1c([H])c([H])c(C([H])=O)c([H])c1O[H].[H]O[H]',
    #   'O=O.CC(O)C(O)C1CNc2[nH]c(=N)nc(O)c2N1.[H]OC(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c([H])c1[H]>>[H]OC(=O)C([H])(N([H])[H])C([H])([H])c1c([H])c([H])c(O[H])c(O[H])c1[H].CC(O)C(O)C1CN=C2NC(=N)N=C(O)C2(O)N1',
    #   '[H]O[H].[H]Oc1c([H])c([H])c(C([H])([H])[H])c([H])c1[H]>>[H]Oc1c([H])c([H])c(C([H])=O)c([H])c1[H]'
    # ]
    rxn_i = 1
    for rxn_id in rxn_lst:
      rxn = pathway.get_reaction(rxn_id)
      data = {}
      data[''] = i
      data['Reaction'] = f'RP{rxn_i}'
      data['Pathway Name'] = pathway.get_id()
      data['Reaction Rule'] = '>>'.join(rxn.get_smiles().split('>>')[::-1])
      data['Rule Score'] = rxn.get_rule_score()
      reactions_fp.write(
        ','.join(
          [str(data.get(c, '')) for c in columns]
        ) + '\n'
      )
      i += 1
      rxn_i += 1

################################################################


