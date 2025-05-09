"""
Source: https://github.com/Coughy1991/Reaction_condition_recommendation/
Using Machine Learning To Predict Suitable Conditions for Organic Reactions, ACS Cent. Sci. 2018, 4, 11, 1465–1476
Authors: Hanyu Gao, Thomas J. Struble, Connor W. Coley, Yuran Wang, William H. Green, Klavs F. Jensen

Modified by: Shivesh Prakash (shivesh.prakash@mail.utoronto.ca)
"""

import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy import stats
from keras import backend as K
from keras.models import model_from_json
from rdkit.Chem import AllChem, DataStructs

contextRecommender_loc = "contextRecommender"


def create_rxn_Morgan2FP_separately(
    rsmi,
    psmi,
    rxnfpsize=16384,
    pfpsize=16384,
    useFeatures=False,
    calculate_rfp=True,
    useChirality=False,
):
    # Similar as the above function but takes smiles separately and returns pfp and rfp separately

    rsmi = rsmi.encode("utf-8")
    psmi = psmi.encode("utf-8")
    try:
        mol = Chem.MolFromSmiles(rsmi)
    except Exception as e:
        print(e)
        return
    try:
        fp_bit = AllChem.GetMorganFingerprintAsBitVect(
            mol=mol,
            radius=2,
            nBits=rxnfpsize,
            useFeatures=useFeatures,
            useChirality=useChirality,
        )
        fp = np.empty(rxnfpsize, dtype="float32")
        DataStructs.ConvertToNumpyArray(fp_bit, fp)
    except Exception as e:
        print("Cannot build reactant fp due to {}".format(e))
        return
    rfp = fp

    try:
        mol = Chem.MolFromSmiles(psmi)
    except Exception as e:
        return
    try:
        fp_bit = AllChem.GetMorganFingerprintAsBitVect(
            mol=mol,
            radius=2,
            nBits=pfpsize,
            useFeatures=useFeatures,
            useChirality=useChirality,
        )
        fp = np.empty(pfpsize, dtype="float32")
        DataStructs.ConvertToNumpyArray(fp_bit, fp)
    except Exception as e:
        print("Cannot build product fp due to {}".format(e))
        return
    pfp = fp
    return [pfp, rfp]


class NeuralNetContextRecommender:
    """Reaction condition predictor based on Nearest Neighbor method"""

    def __init__(self, max_contexts=10, singleSlvt=True, with_smiles=True):
        """
        :param singleSlvt:
        :param with_smiles:
        """
        self.nnModel = None
        self.c1_dict = None
        self.s1_dict = None
        self.s2_dict = None
        self.r1_dict = None
        self.r2_dict = None
        self.num_cond = 1
        self.singleSlvt = singleSlvt
        self.with_smiles = with_smiles
        self.max_total_context = max_contexts
        self.max_context = 2
        self.fp_size = 2048

    def load(self, model_path="", info_path="", weights_path=""):
        # for the neural net model, info path points to the encoders
        self.load_nn_model(model_path, info_path, weights_path)

        print("Nerual network context recommender has been loaded.")

    def load_nn_model(self, model_path="", info_path="", weights_path=""):

        if not model_path:
            print(
                "Cannot load neural net context recommender without a specific path to the model. Exiting..."
            )
        if not info_path:
            print(
                "Cannot load nerual net context recommender without a specific path to the model info. Exiting..."
            )

        ###load model##############
        # load json and create model
        K.set_learning_phase(0)
        json_file = open(model_path, "r")
        loaded_model_json = json_file.read()
        json_file.close()

        self.nnModel = model_from_json(loaded_model_json)
        # load weights into new model
        self.nnModel.load_weights(weights_path)
        # get fp_size based on the model
        self.fp_size = self.nnModel.input_shape[0][1]
        r1_dict_file = info_path + "r1_dict.pickle"
        r2_dict_file = info_path + "r2_dict.pickle"
        s1_dict_file = info_path + "s1_dict.pickle"
        s2_dict_file = info_path + "s2_dict.pickle"
        c1_dict_file = info_path + "c1_dict.pickle"

        with open(r1_dict_file, "rb") as R1_DICT_F:
            self.r1_dict = pickle.load(R1_DICT_F)

        with open(r2_dict_file, "rb") as R2_DICT_F:
            self.r2_dict = pickle.load(R2_DICT_F)

        with open(s1_dict_file, "rb") as S1_DICT_F:
            self.s1_dict = pickle.load(S1_DICT_F)

        with open(s2_dict_file, "rb") as S2_DICT_F:
            self.s2_dict = pickle.load(S2_DICT_F)

        with open(c1_dict_file, "rb") as C1_DICT_F:
            self.c1_dict = pickle.load(C1_DICT_F)

        self.inverse_c1_dict = {v: k for k, v in self.c1_dict.items()}
        self.inverse_s1_dict = {v: k for k, v in self.s1_dict.items()}
        self.inverse_s2_dict = {v: k for k, v in self.s2_dict.items()}
        self.inverse_r1_dict = {v: k for k, v in self.r1_dict.items()}
        self.inverse_r2_dict = {v: k for k, v in self.r2_dict.items()}

        self.c1_dim = self.nnModel.input_shape[2][1]
        self.r1_dim = self.nnModel.input_shape[3][1]
        self.r2_dim = self.nnModel.input_shape[4][1]
        self.s1_dim = self.nnModel.input_shape[5][1]
        self.s2_dim = self.nnModel.input_shape[6][1]
        fp_transform_layer = self.nnModel.get_layer("fp_transform1")
        c1_input_layer = self.nnModel.get_layer("input_c1")
        s1_input_layer = self.nnModel.get_layer("input_s1")
        s2_input_layer = self.nnModel.get_layer("input_s2")
        r1_input_layer = self.nnModel.get_layer("input_r1")
        r2_input_layer = self.nnModel.get_layer("input_r2")
        c1_output = self.nnModel.get_layer("c1")
        s1_output = self.nnModel.get_layer("s1")
        s2_output = self.nnModel.get_layer("s2")
        r1_output = self.nnModel.get_layer("r1")
        r2_output = self.nnModel.get_layer("r2")
        T_output = self.nnModel.get_layer("T")

        self.fp_func = K.function(self.nnModel.inputs, [fp_transform_layer.output])
        self.c1_func = K.function([fp_transform_layer.output], [c1_output.output])
        self.s1_func = K.function(
            [fp_transform_layer.output, c1_input_layer.output], [s1_output.output]
        )
        self.s2_func = K.function(
            [fp_transform_layer.output, c1_input_layer.output, s1_input_layer.output],
            [s2_output.output],
        )
        self.r1_func = K.function(
            [
                fp_transform_layer.output,
                c1_input_layer.output,
                s1_input_layer.output,
                s2_input_layer.output,
            ],
            [r1_output.output],
        )
        self.r2_func = K.function(
            [
                fp_transform_layer.output,
                c1_input_layer.output,
                s1_input_layer.output,
                s2_input_layer.output,
                r1_input_layer.output,
            ],
            [r2_output.output],
        )
        self.T_func = K.function(
            [
                fp_transform_layer.output,
                c1_input_layer.output,
                s1_input_layer.output,
                s2_input_layer.output,
                r1_input_layer.output,
                r2_input_layer.output,
            ],
            [T_output.output],
        )

    def load_predictor(self, userInput):
        """Loads the predictor based on user input"""
        self.num_cond = userInput["num_cond"]
        self.dist_limit = userInput["dist_limit"]
        self.singleSlvt = userInput["first_solvent_only"]
        self.with_smiles = userInput["with_smiles_only"]
        self.max_total_context = userInput["max_total_context"]
        self.max_int = userInput["max_int"]
        self.max_context = userInput["max_context"]

    def get_n_conditions(
        self, rxn, n=10, singleSlvt=True, with_smiles=True, return_scores=False
    ):
        """
        Reaction condition recommendations for a rxn (SMILES) from top n NN
        Returns the top n parseable conditions.
        """
        # print('started neuralnet')
        self.singleSlvt = singleSlvt
        self.with_smiles = with_smiles
        # print('rxn to recommend context for : {}'.format(rxn), contextRecommender_loc)
        try:
            rsmi = rxn.split(">>")[0]
            psmi = rxn.split(">>")[1]
            rct_mol = Chem.MolFromSmiles(rsmi)
            prd_mol = Chem.MolFromSmiles(psmi)
            [
                atom.ClearProp("molAtomMapNumber")
                for atom in rct_mol.GetAtoms()
                if atom.HasProp("molAtomMapNumber")
            ]
            [
                atom.ClearProp("molAtomMapNumber")
                for atom in prd_mol.GetAtoms()
                if atom.HasProp("molAtomMapNumber")
            ]
            rsmi = Chem.MolToSmiles(rct_mol, isomericSmiles=True)
            psmi = Chem.MolToSmiles(prd_mol, isomericSmiles=True)
            [pfp, rfp] = create_rxn_Morgan2FP_separately(
                rsmi,
                psmi,
                rxnfpsize=self.fp_size,
                pfpsize=self.fp_size,
                useFeatures=False,
                calculate_rfp=True,
                useChirality=True,
            )
            pfp = pfp.reshape(1, self.fp_size)
            rfp = rfp.reshape(1, self.fp_size)
            rxnfp = pfp - rfp
            c1_input = []
            r1_input = []
            r2_input = []
            s1_input = []
            s2_input = []
            inputs = [pfp, rxnfp, c1_input, r1_input, r2_input, s1_input, s2_input]

            (top_combos, top_combo_scores) = self.predict_top_combos(inputs=inputs)

            if return_scores:
                return (top_combos[:n], top_combo_scores[:n])
            else:
                return top_combos[:n]

        except Exception as e:
            print(
                "Failed for reaction {} because {}. Returning None.".format(rxn, e),
                contextRecommender_loc,
            )
            return [[]]

    def path_condition(self, n, path):
        """Reaction condition recommendation for a reaction path with multiple reactions
        path: a list of reaction SMILES for each step
        return: a list of reaction context with n options for each step
        """
        rsmi_list = []
        psmi_list = []
        contexts = []
        for rxn in path:
            try:
                rsmi = rxn.split(">>")[0]
                psmi = rxn.split(">>")[1]

                rct_mol = Chem.MolFromSmiles(rsmi)
                prd_mol = Chem.MolFromSmiles(psmi)
                [
                    atom.ClearProp("molAtomMapNumber")
                    for atom in rct_mol.GetAtoms()
                    if atom.HasProp("molAtomMapNumber")
                ]
                [
                    atom.ClearProp("molAtomMapNumber")
                    for atom in prd_mol.GetAtoms()
                    if atom.HasProp("molAtomMapNumber")
                ]
                rsmi = Chem.MolToSmiles(rct_mol, isomericSmiles=True)
                psmi = Chem.MolToSmiles(prd_mol, isomericSmiles=True)
                [pfp, rfp] = create_rxn_Morgan2FP_separately(
                    rsmi,
                    psmi,
                    rxnfpsize=self.fp_size,
                    pfpsize=self.fp_size,
                    useFeatures=False,
                    calculate_rfp=True,
                    useChirality=True,
                )
                pfp = pfp.reshape(1, self.fp_size)
                rfp = rfp.reshape(1, self.fp_size)
                rxnfp = pfp - rfp
                c1_input = []
                r1_input = []
                r2_input = []
                s1_input = []
                s2_input = []
                inputs = [pfp, rxnfp, c1_input, r1_input, r2_input, s1_input, s2_input]
                top_combos = self.predict_top_combos(
                    inputs=inputs,
                    c1_rank_thres=1,
                    s1_rank_thres=3,
                    s2_rank_thres=1,
                    r1_rank_thres=4,
                    r2_rank_thres=1,
                )
                contexts.append(top_combos[:n])
            except Exception as e:
                print(
                    "Failed for reaction {} because {}. Returning None.".format(rxn, e)
                )
        return contexts

    def predict_top_combos(
        self,
        inputs,
        return_categories_only=False,
        c1_rank_thres=2,
        s1_rank_thres=3,
        s2_rank_thres=1,
        r1_rank_thres=3,
        r2_rank_thres=1,
    ):
        # this function predicts the top combos based on rank thresholds for
        # individual elements
        context_combos = []
        context_combo_scores = []
        num_combos = (
            c1_rank_thres
            * s1_rank_thres
            * s2_rank_thres
            * r1_rank_thres
            * r2_rank_thres
        )
        [
            pfp,
            rxnfp,
            c1_input_user,
            r1_input_user,
            r2_input_user,
            s1_input_user,
            s2_input_user,
        ] = inputs

        # set s2 to be none ## specifically for the single sovlent case
        # s2_input_user = np.zeros(s2_dim,dtype = 'float32').reshape(1,s2_dim)
        # s2_input_user[0] = 1
        c1_input_dum = np.zeros(self.c1_dim, dtype="float32").reshape(1, self.c1_dim)
        r1_input_dum = np.zeros(self.r1_dim, dtype="float32").reshape(1, self.r1_dim)
        r2_input_dum = np.zeros(self.r2_dim, dtype="float32").reshape(1, self.r2_dim)
        s1_input_dum = np.zeros(self.s1_dim, dtype="float32").reshape(1, self.s1_dim)
        s2_input_dum = np.zeros(self.s2_dim, dtype="float32").reshape(1, self.s2_dim)
        model_inputs = [
            pfp,
            rxnfp,
            c1_input_dum,
            r1_input_dum,
            r2_input_dum,
            s1_input_dum,
            s2_input_dum,
        ]
        fp_trans = self.fp_func(model_inputs)
        if c1_input_user == []:
            c1_inputs = fp_trans
            c1_pred = self.c1_func(c1_inputs)
            c1_cdts = c1_pred[0][0].argsort()[-c1_rank_thres:][::-1]
        else:
            c1_cdts = np.nonzero(c1_input_user)[0]
        # find the name of catalyst
        for c1_cdt in c1_cdts:
            c1_name = self.c1_dict[c1_cdt]
            c1_input = np.zeros([1, self.c1_dim])
            c1_input[0, c1_cdt] = 1
            if c1_input_user == []:
                c1_sc = c1_pred[0][0][c1_cdt]
            else:
                c1_sc = 1
            if s1_input_user == []:
                s1_inputs = [fp_trans[0], c1_input]
                s1_pred = self.s1_func(s1_inputs)
                s1_cdts = s1_pred[0][0].argsort()[-s1_rank_thres:][::-1]
            else:
                s1_cdts = np.nonzero(s1_input_user)[0]
            for s1_cdt in s1_cdts:
                s1_name = self.s1_dict[s1_cdt]
                s1_input = np.zeros([1, self.s1_dim])
                s1_input[0, s1_cdt] = 1
                if s1_input_user == []:
                    s1_sc = s1_pred[0][0][s1_cdt]
                else:
                    s1_sc = 1
                if s2_input_user == []:
                    s2_inputs = [fp_trans[0], c1_input, s1_input]
                    s2_pred = self.s2_func(s2_inputs)
                    s2_cdts = s2_pred[0][0].argsort()[-s2_rank_thres:][::-1]
                else:
                    s2_cdts = np.nonzero(s2_input_user)[0]
                for s2_cdt in s2_cdts:
                    s2_name = self.s2_dict[s2_cdt]
                    s2_input = np.zeros([1, self.s2_dim])
                    s2_input[0, s2_cdt] = 1
                    if s2_input_user == []:
                        s2_sc = s2_pred[0][0][s2_cdt]
                    else:
                        s2_sc = 1
                    if r1_input_user == []:
                        r1_inputs = [fp_trans[0], c1_input, s1_input, s2_input]
                        r1_pred = self.r1_func(r1_inputs)
                        r1_cdts = r1_pred[0][0].argsort()[-r1_rank_thres:][::-1]
                    else:
                        r1_cdts = np.nonzero(r1_input_user)[0]
                    for r1_cdt in r1_cdts:
                        r1_name = self.r1_dict[r1_cdt]
                        r1_input = np.zeros([1, self.r1_dim])
                        r1_input[0, r1_cdt] = 1
                        if r1_input_user == []:
                            r1_sc = r1_pred[0][0][r1_cdt]
                        else:
                            r1_sc = 1
                        if r2_input_user == []:
                            r2_inputs = [
                                fp_trans[0],
                                c1_input,
                                s1_input,
                                s2_input,
                                r1_input,
                            ]
                            r2_pred = self.r2_func(r2_inputs)
                            r2_cdts = r2_pred[0][0].argsort()[-r2_rank_thres:][::-1]
                        else:
                            r2_cdts = np.nonzero(r2_input_user)[0]
                        for r2_cdt in r2_cdts:
                            r2_name = self.r2_dict[r2_cdt]
                            r2_input = np.zeros([1, self.r2_dim])
                            r2_input[0, r2_cdt] = 1
                            if r2_input_user == []:
                                r2_sc = r2_pred[0][0][r2_cdt]
                            else:
                                r2_sc = 1
                            T_inputs = [
                                fp_trans[0],
                                c1_input,
                                s1_input,
                                s2_input,
                                r1_input,
                                r2_input,
                            ]
                            T_pred = self.T_func(T_inputs)
                            # print(c1_name,s1_name,s2_name,r1_name,r2_name)
                            cat_name = [c1_name]
                            if r2_name == "":
                                rgt_name = [r1_name]
                            else:
                                rgt_name = [r1_name, r2_name]
                            if s2_name == "":
                                slv_name = [s1_name]
                            else:
                                slv_name = [s1_name, s2_name]
                            # if self.with_smiles:
                            #     rgt_name = [rgt for rgt in rgt_name if 'Reaxys' not in rgt]
                            #     slv_name = [slv for slv in slv_name if 'Reaxys' not in slv]
                            #     cat_name = [cat for cat in cat_name if 'Reaxys' not in cat]
                            ##for testing purpose only, output order as training
                            if return_categories_only:
                                context_combos.append(
                                    [
                                        c1_cdt,
                                        s1_cdt,
                                        s2_cdt,
                                        r1_cdt,
                                        r2_cdt,
                                        T_pred[0][0][0],
                                    ]
                                )
                            ## esle ouptupt format compatible with the overall framework
                            else:
                                context_combos.append(
                                    [
                                        float(T_pred[0][0][0]),
                                        ".".join(slv_name),
                                        ".".join(rgt_name),
                                        ".".join(cat_name),
                                        np.nan,
                                        np.nan,
                                    ]
                                )

                            context_combo_scores.append(
                                c1_sc * s1_sc * s2_sc * r1_sc * r2_sc
                            )
        context_ranks = list(num_combos + 1 - stats.rankdata(context_combo_scores))

        context_combos = [
            context_combos[context_ranks.index(i + 1)] for i in range(num_combos)
        ]
        context_combo_scores = [
            context_combo_scores[context_ranks.index(i + 1)] for i in range(num_combos)
        ]

        return (context_combos, context_combo_scores)

    def category_to_name(self, chem_type, category):
        if chem_type == "c1":
            return self.c1_dict[category]
        elif chem_type == "s1":
            return self.s1_dict[category]
        elif chem_type == "s2":
            return self.s2_dict[category]
        elif chem_type == "r1":
            return self.r1_dict[category]
        elif chem_type == "r2":
            return self.r2_dict[category]

    def name_to_category(self, chem_type, name):
        try:
            if chem_type == "c1":
                return self.inverse_c1_dict[name]
            elif chem_type == "s1":
                return self.inverse_s1_dict[name]
            elif chem_type == "s2":
                return self.inverse_s2_dict[name]
            elif chem_type == "r1":
                return self.inverse_r1_dict[name]
            elif chem_type == "r2":
                return self.inverse_r2_dict[name]
        except:
            print("name {} not found!".format(name))


def pred_temperature(rxn):
    """
    Predict optimal reaction temperature using neural network recommendations.

    Calculates a weighted average temperature from top 10 model predictions,
    with weights based on prediction confidence scores.
    """
    cont = NeuralNetContextRecommender()
    cont.load_nn_model(
        model_path="NeuralNet_Cont_Model/model.json",
        info_path="NeuralNet_Cont_Model/",
        weights_path="NeuralNet_Cont_Model/weights.h5",
    )
    items, scores = cont.get_n_conditions(
        rxn, 10, with_smiles=False, return_scores=True
    )
    temp = 0
    total = sum(scores)
    for i in range(10):
        temp += items[i][0] * scores[i] / total
    return temp


def pred_solvent_score(rxn, solvent_path="data/toxicity.csv"):
    """
    Calculate environmental toxicity score for recommended solvents.

    Computes weighted average toxicity score using top 10 solvent recommendations,
    combining both primary and secondary solvent scores from a toxicity database.
    """
    cont = NeuralNetContextRecommender()
    cont.load_nn_model(
        model_path="NeuralNet_Cont_Model/model.json",
        info_path="NeuralNet_Cont_Model/",
        weights_path="NeuralNet_Cont_Model/weights.h5",
    )
    items, scores = cont.get_n_conditions(
        rxn, 10, with_smiles=False, return_scores=True
    )
    solvent_score = 0
    total = sum(scores)
    df = pd.read_csv(solvent_path)
    for i in range(10):
        smile = items[i][1]
        row = df[df["smiles"] == smile]
        if not row.empty:
            solvent_score += int(row["score"]) * scores[i] / total
        smile = items[i][2]
        row = df[df["smiles"] == smile]
        if not row.empty:
            solvent_score += int(row["score"]) * scores[i] / total
    return solvent_score / 2


def get_solvent_reagent(rxn):
    """
    Retrieve top recommended solvent and reagent combination for a reaction.

    Returns the most confidently predicted solvent-reagent pair from the
    neural network model's top recommendation.
    """
    try:
        cont = NeuralNetContextRecommender()
        cont.load_nn_model(
            model_path="NeuralNet_Cont_Model/model.json",
            info_path="NeuralNet_Cont_Model/",
            weights_path="NeuralNet_Cont_Model/weights.h5",
        )
        items, scores = cont.get_n_conditions(
            rxn, 10, with_smiles=False, return_scores=True
        )
        return items[0][1], items[0][2]
    except:
        return None, None
