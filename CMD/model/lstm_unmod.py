from __future__ import absolute_import
from __future__ import print_function
import numpy as np
np.random.seed(1337) # for reproducibility

import time
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
from scipy import stats
from scipy.stats.stats import pearsonr

from keras.preprocessing import sequence
from keras.optimizers import SGD, RMSprop, Adagrad
from keras.utils import np_utils
from keras.models import Sequential
from keras.models import load_model
from keras.layers.core import Dense, Dropout, Activation, TimeDistributedDense, Masking
from keras.layers.embeddings import Embedding
from keras.layers.recurrent import LSTM, GRU, SimpleRNN
from keras.layers.wrappers import TimeDistributed
from keras.layers.wrappers import Bidirectional
from keras.initializations import normal, identity

from .featurize import Ion2Vector

import fnmatch
import os

class IonLSTM:
    # conf is ETD_config or HCD_Config in fragmentation_config.py
    def __init__(self, conf):
        self.time_step = conf.time_step # max peptide_length is 20, so there are at most 19 ion sites, LSTM time_step should be peptide_length - 1
        self.max_charge = conf.max_ion_charge
        self.ion_types = conf.ion_types
        self.from_CTerm = conf.from_CTerm
        self.batch_size = 1024
        self.layer_size = 256
        self.epoches = 100
        self.config = conf
        self.step_epoches = 20 # for TrainWithStep() only
        self.initial_epoch = 0
        self.layer_count = 0
        self.retrain = True
        self.goontrain = False
        self.model = None
        self.debug = False
        if self.debug:
            self.layer_size = 256
            self.epoches = 100
            self.step_epoches = 20 # for TrainWithStep() only

    # compare the predicted y with real y, for LSTM only
    def ComparePredReal_RNN(self, predictions, real, peplens):
        print("    pred shape", predictions.shape)
        ypred_seq = np.reshape(predictions, (predictions.shape[0], predictions.shape[1] * predictions.shape[2]), order='C')
        ytest_seq = np.reshape(real, (real.shape[0], real.shape[1] * real.shape[2]), order='C')
        sims = []
        # lenset = set(peplens)
        # peplens = np.array(peplens)
        # for plen in lenset:
            # ypred_tmp = ypred_seq[peplens == plen, :(plen-1)*predictions.shape[2]]
            # ytest_tmp = ytest_seq[peplens == plen, :(plen-1)*predictions.shape[2]]
            # for x, y in np.nditer([ypred_tmp, ytest_tmp]):
                # sim = pearsonr(x, y)[0]
                # if not np.isnan(sim):
                    # sims.append(sim)
        for i in range(len(predictions)):
            # sim = pearsonr(np.reshape(predictions[i,:(peplens[i]-1),:],-1), np.reshape(real[i,:(peplens[i]-1),:],-1))[0]
            sim = pearsonr(ypred_seq[i][:(peplens[i]-1) * predictions.shape[2]], ytest_seq[i][:(peplens[i]-1) * predictions.shape[2]])[0]
            sims.append(sim)
        sims_nan = np.array(sims)
        sims = sims_nan[np.isnan(sims_nan) == False]
        med = np.median(sims)
        mad = np.median(np.abs(sims - med))
        avg = np.mean(sims)
        std = np.std(sims)
        out_median = "    Median sim = %.3f, MAD sim = %.3f" %(med, mad)
        out_mean = "    Mean sim = %.3f, STD sim = %.3f" %(avg, std)
        print(out_median)
        print(out_mean)
        return (sims_nan,(med,mad,avg,std))
    
    def ComparePredReal_notRNN(self, predictions, real, peplens):
        print("    pred shape", predictions.shape)
        sims = []
        start = 0
        for i in range(len(peplens)):
            sim = pearsonr(predictions[start:(start + (peplens[i]-1) * predictions.shape[2])], real[start:(start + (peplens[i]-1) * predictions.shape[2])])[0]
            if not np.isnan(sim):
                sims.append(sim)
            else:
                sims.append(0)
            start += (peplens[i]-1) * predictions.shape[2]
        
        med = np.median(sims)
        mad = np.median(np.abs(sims - med))
        avg = np.mean(sims)
        std = np.std(sims)
        out_median = "    Median sim = %.3f, MAD sim = %.3f" %(med, mad)
        out_mean = "    Mean sim = %.3f, STD sim = %.3f\n" %(avg, std)
        print(out_median)
        print(out_mean)
        return (sims,(med,mad,avg,std))

    # data
    def DirectIon2Vector(self, filenames):
        ion2vec = Ion2Vector(modini = self.config.modini, ion_types = self.ion_types, cleave_site_num = self.time_step, max_charge = self.max_charge, prev = 1, next = 1, from_CTerm = self.from_CTerm)
        X = []
        Y = []
        peplen = []
        for filename in filenames:
            print(filename)
            _X,_Y,_peplen = ion2vec.Featurize(filename, RNN = True)
            print (np.array(_X).shape, np.array(_Y).shape)
            X.extend(_X)
            Y.extend(_Y)
            peplen.extend(_peplen)
            
        X = np.array(X)
        Y = np.array(Y)
        peplen = np.array(peplen)
        print (X.shape, Y.shape)
        return X,Y,peplen
    
    # model
    def BuildModel(self, input_shape, nlayers = 1):
        print('BuildModel ... ')
        self.layer_count = nlayers
        self.model = Sequential()
        self.model.add(Masking(mask_value = 0, input_shape = input_shape))
        
        for i in range(nlayers):
            self.model.add(Bidirectional(LSTM(self.layer_size, return_sequences = True, dropout_W = 0.2, dropout_U = 0.2)))
            # self.model.add(LSTM(self.layer_size, return_sequences = True))
            self.model.add(Activation('relu'))
            self.model.add(Dropout(0.3))
            
            # self.model.add(TimeDistributed(Dense(self.layer_size)))
            # self.model.add(Activation('relu'))
            # self.model.add(Dropout(0.3))
        
        # self.model.add(Bidirectional(LSTM(self.layer_size, return_sequences = True, dropout_W = 0.2, dropout_U = 0.2)))
        # self.model.add(Activation('relu'))
        # self.model.add(Dropout(0.3))
        # self.layer_count += 1
        
        # self.model.add(TimeDistributed(Dense(self.layer_size)))
        # self.model.add(Activation('relu'))
        # self.model.add(Dropout(0.3))
        # self.layer_count += 1
        
        #the output layer
        self.model.add(TimeDistributed(Dense(self.max_charge * len(self.ion_types))))
        self.model.add(Activation('relu'))
        
    def InitModel(self, optimizer = 'adam', init_weights_file = None):
        if optimizer == 'sgd':
            optimizer = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        
        # in keras/objectives.py
        # mean_squared_error / mse
        # root_mean_squared_error / rmse
        # mean_absolute_error / mae
        self.model.compile(loss='mae', optimizer=optimizer)
        if not init_weights_file is None:
            self.model.load_weights(init_weights_file)
    
    def TrainModel(self, X_train, y_train, valid = None):
        self.model.fit(X_train, y_train, batch_size=self.batch_size, nb_epoch=self.epoches+self.initial_epoch, validation_split=0.0, initial_epoch=self.initial_epoch)

    def SaveWeights(self, weights_file):
        self.model.save_weights(weights_file)
    # model
    def TestModel(self, X_test, y_test, peplen_test, X_train = None, y_train = None, peplen_train = None):
        toplot = []
        if X_train is not None:
            print("Training set sim:")
            predictions = self.ModelPredict(X_train)
            compRNN = self.ComparePredReal_RNN(predictions, y_train, peplen_train)
            toplot.append(compRNN[0])
            ResTrain = compRNN[1]
        else:
            ResTrain = (0,0,0,0)
        print("Test set sim:")
        predictions = self.ModelPredict(X_test)
        compRNN = self.ComparePredReal_RNN(predictions, y_test, peplen_test)
        toplot.append(compRNN[0])
        ResTest = compRNN[1]
        plt.figure()
        plt.boxplot(toplot,sym='')
        plt.savefig("./models/lstm_layer=%d-epoch=%d-ch=%d-type=%s.png" %(self.layer_count, self.epoches+self.initial_epoch, self.max_charge, "".join(self.ion_types)))
        return (ResTrain,ResTest)

    # model
    def Predict(self, X):
        predictions = self.model.predict(X, batch_size = self.batch_size)
        predictions[predictions < 0] = 0
        predictions[predictions > 1] = 1
        return predictions

    # model
    def SaveModel(self, filepath = None):
        if filepath is None:
            self.model.save("./models/lstm_layer=%d-epoch=%d-ch=%d-type=%s.lstm" %(self.layer_count, self.epoches+self.initial_epoch, self.max_charge, "".join(self.ion_types)))
        else:
            self.model.save(filepath)
    
    # model
    def LoadModel(self, modelfile):
        self.model = load_model(modelfile)

    # data
    def LoadDataFolder(self, dataset_folder):
        print("Loading data from folder...")
        filenames = []
        for input_file in os.listdir(dataset_folder):
            if fnmatch.fnmatch(input_file, "*.plabel"):
                filenames.append(os.path.join(dataset_folder, input_file))
        return self.DirectIon2Vector(filenames)

    # data
    def LoadDataFiles(self, filenames):
        print("Loading data from files...")
        return self.DirectIon2Vector(filenames)
        
    def SetDataFiles(train_files, test_files):
        self.train_filenames = train_files
        self.test_filenames = test_files
    
    def SetDataFolder(train_folder, test_folder):
        self.traindata_folder = train_folder
        self.testdata_folder = test_folder
    
    def Run(self):
        # main
        self.traindata_folder = "./Mann-MouseBrain-NNeu-2015-QEHF-27"
        self.testdata_folder = "../data/Mann-HeLa-MCP-2011-QE-30"
        
        # train_filenames = ["../data/Mann-HeLa-MCP-2011-Velos-40/20100825_Velos2_AnMi_QC_wt_HCD_iso4_swG_HCDFT.plabel",
                           # "./Mann-HeLa-MCP-2011-Velos-40/20100825_Velos2_AnMi_QC_wt_HCD_iso4_swG_2_HCDFT.plabel"]
        # test_filenames = ["./Mann-Synthetic-NBT-2013-Velos-40/5_HCDFT.plabel"]
        
        self.train_filenames = ["../data/Mann-MouseBrain-NNeu-2015-QEHF-27/20141202_QEp8_KiSh_SA_Cerebellum_P05_Singleshot1_HCDFT.plabel"]
        # test_filenames = ["./Mann-HeLa-MCP-2011-QE-30/20110405_AnMi_HeLa_matching_HCDFT.plabel"]
        # test_filenames = ["./Mann-HeLa-MCP-2011-QE-30/20110405_AnMi_HeLa_matching_HCDFT.plabel.ch2.txt"]
        if self.retrain or self.goontrain:
            if self.debug:
                X_train, y_train, peplen_train = self.LoadDataFiles(self.train_filenames)
            else:
                X_train, y_train, peplen_train = self.LoadDataFolder(self.traindata_folder)
            print(len(X_train), 'train sequences')
        else:
            X_train = None
            y_train = None
            peplen_train = None

        X_test, y_test, peplen_test = self.LoadDataFolder(self.testdata_folder)
        # X_test, y_test, peplen_test = self.LoadDataFiles(self.test_filenames)
        print(len(X_test), 'test sequences')
        
        self.model = None
        if self.retrain:
            self.BuildModel(input_shape = (X_train.shape[1], X_train.shape[2]),nlayers=3)
            self.InitModel(optimizer = 'adam')
            self.TrainModel(X_train, y_train)
            self.SaveModel()
        elif self.goontrain:
            print("Resume training ...")
            self.LoadModel(modelfile = self.modelfile)
            self.model.fit(X_train, y_train, batch_size=self.batch_size, nb_epoch=self.epoches+self.initial_epoch, validation_split=0.0, initial_epoch=self.initial_epoch)
        else:
            self.LoadModel(modelfile = self.modelfile)
        self.TestModel(X_test, y_test, peplen_test, X_train, y_train, peplen_train)
    
    def TrainWithStep(self):
        self.ccfile = "pearson_cc.txt"
        outcc = open(self.ccfile,"w")
        outcc.close()
        traindata_folder = "./Mann-MouseBrain-NNeu-2015-QEHF-27"
        testdata_folder = "./Mann-HeLa-MCP-2011-QE-30"
        train_filenames = ["./Mann-MouseBrain-NNeu-2015-QEHF-27/20141202_QEp8_KiSh_SA_Cerebellum_P05_Singleshot1_HCDFT.plabel"]
        test_filenames = ["./Mann-HeLa-MCP-2011-QE-30/20110405_AnMi_HeLa_matching_HCDFT.plabel"]
        # test_filenames = ["./Mann-HeLa-MCP-2011-QE-30/20110405_AnMi_HeLa_matching_HCDFT.plabel.ch2.txt"]
        if self.debug:
            X_train, y_train, peplen_train = self.LoadDataFiles(train_filenames)
            X_test, y_test, peplen_test = self.LoadDataFiles(test_filenames)
        else:
            X_train, y_train, peplen_train = self.LoadDataFolder(traindata_folder)
            X_test, y_test, peplen_test = self.LoadDataFolder(testdata_folder)
        print(len(X_train), 'train sequences')
        print(len(X_test), 'test sequences')
        
        self.final_epoches = self.epoches
        self.epoches = self.step_epoches
        self.initial_epoch = 0
        
        weights_file = 'lstm_adam.weights'
        
        self.BuildModel(input_shape = (X_train.shape[1], X_train.shape[2]))
        self.InitModel(optimizer = 'adam')
        self.TrainModel(X_train, y_train)
        self.SaveWeights(weights_file)
        self.SaveModel()
        self.outcc = open(self.ccfile,"a")
        self.outcc.write("#epoch | med CC (tr) | mean CC (tr) | med CC (ts) | mean CC (ts)\n-----|----------:|----------:|----------:|----------:\n")
        (ResTrain, ResTest) = self.TestModel(X_test, y_test, peplen_test, X_train, y_train, peplen_train)
        self.outcc.write("%d | %.3f | %.3f | %.3f | %.3f\n" %(self.epoches + self.initial_epoch, ResTrain[0], ResTrain[2], ResTest[0], ResTest[2]))
        self.outcc.close()
        
        self.initial_epoch = self.epoches
        self.BuildModel(input_shape = (X_train.shape[1], X_train.shape[2]))
        self.InitModel(optimizer = 'sgd', init_weights_file = weights_file)
        self.TrainModel(X_train, y_train)
        self.SaveModel()
        self.outcc = open(self.ccfile,"a")
        (ResTrain, ResTest) = self.TestModel(X_test, y_test, peplen_test, X_train, y_train, peplen_train)
        self.outcc.write("%d | %.3f | %.3f | %.3f | %.3f\n" %(self.epoches + self.initial_epoch, ResTrain[0], ResTrain[2], ResTest[0], ResTest[2]))
        self.outcc.close()
        
        for self.initial_epoch in range(self.epoches, self.final_epoches, self.step_epoches):
            self.TrainModel(X_train, y_train)
            print("epoches = %d" %(self.epoches+self.initial_epoch))
            self.SaveModel()
            self.outcc = open(self.ccfile,"a")
            (ResTrain, ResTest) = self.TestModel(X_test, y_test, peplen_test, X_train, y_train, peplen_train)
            self.outcc.write("%d | %.3f | %.3f | %.3f | %.3f \n" %(self.epoches + self.initial_epoch, ResTrain[0], ResTrain[2], ResTest[0], ResTest[2]))
            self.outcc.close()
            

class IonLSTM_seq_level:
    # conf is ETD_config or HCD_Config in fragmentation_config.py
    def __init__(self, conf):
        self.time_step = conf.time_step # max peptide_length is 20, so there are at most 19 ion sites, LSTM time_step should be peptide_length - 1
        self.max_charge = conf.max_ion_charge
        self.ion_types = conf.ion_types
        self.from_CTerm = conf.from_CTerm
        self.batch_size = 1024
        self.layer_size = 256
        self.epoches = 100
        self.config = conf
        self.step_epoches = 20 # for TrainWithStep() only
        self.initial_epoch = 0
        self.layer_count = 0
        self.retrain = True
        self.goontrain = False
        self.model = None
        self.debug = False
        if self.debug:
            self.layer_size = 256
            self.epoches = 100
            self.step_epoches = 20 # for TrainWithStep() only

    # compare the predicted y with real y, for LSTM only
    def ComparePredReal_RNN(self, predictions, real, peplens):
        print("    pred shape", predictions.shape)
        ypred_seq = np.reshape(predictions, (predictions.shape[0], predictions.shape[1] * predictions.shape[2]), order='C')
        ytest_seq = np.reshape(real, (real.shape[0], real.shape[1] * real.shape[2]), order='C')
        sims = []
        # lenset = set(peplens)
        # peplens = np.array(peplens)
        # for plen in lenset:
            # ypred_tmp = ypred_seq[peplens == plen, :(plen-1)*predictions.shape[2]]
            # ytest_tmp = ytest_seq[peplens == plen, :(plen-1)*predictions.shape[2]]
            # for x, y in np.nditer([ypred_tmp, ytest_tmp]):
                # sim = pearsonr(x, y)[0]
                # if not np.isnan(sim):
                    # sims.append(sim)
        for i in range(len(predictions)):
            # sim = pearsonr(np.reshape(predictions[i,:(peplens[i]-1),:],-1), np.reshape(real[i,:(peplens[i]-1),:],-1))[0]
            sim = pearsonr(ypred_seq[i][:(peplens[i]-1) * predictions.shape[2]], ytest_seq[i][:(peplens[i]-1) * predictions.shape[2]])[0]
            sims.append(sim)
        sims_nan = np.array(sims)
        sims = sims_nan[np.isnan(sims_nan) == False]
        med = np.median(sims)
        mad = np.median(np.abs(sims - med))
        avg = np.mean(sims)
        std = np.std(sims)
        out_median = "    Median sim = %.3f, MAD sim = %.3f" %(med, mad)
        out_mean = "    Mean sim = %.3f, STD sim = %.3f" %(avg, std)
        print(out_median)
        print(out_mean)
        return (sims_nan,(med,mad,avg,std))
    
    def ComparePredReal_notRNN(self, predictions, real, peplens):
        print("    pred shape", predictions.shape)
        sims = []
        start = 0
        for i in range(len(peplens)):
            sim = pearsonr(predictions[start:(start + (peplens[i]-1) * predictions.shape[2])], real[start:(start + (peplens[i]-1) * predictions.shape[2])])[0]
            if not np.isnan(sim):
                sims.append(sim)
            else:
                sims.append(0)
            start += (peplens[i]-1) * predictions.shape[2]
        
        med = np.median(sims)
        mad = np.median(np.abs(sims - med))
        avg = np.mean(sims)
        std = np.std(sims)
        out_median = "    Median sim = %.3f, MAD sim = %.3f" %(med, mad)
        out_mean = "    Mean sim = %.3f, STD sim = %.3f\n" %(avg, std)
        print(out_median)
        print(out_mean)
        return (sims,(med,mad,avg,std))

    # data
    def DirectIon2Vector(self, filenames):
        ion2vec = Ion2Vector(modini = self.config.modini, ion_types = self.ion_types, cleave_site_num = self.time_step, max_charge = self.max_charge, prev = 1, next = 1, from_CTerm = self.from_CTerm)
        X = []
        Y = []
        peplen = []
        for filename in filenames:
            print(filename)
            _X,_Y,_peplen = ion2vec.Featurize_ignore_mod(filename, RNN = True)
            print (np.array(_X).shape, np.array(_Y).shape)
            X.extend(_X)
            Y.extend(_Y)
            peplen.extend(_peplen)
            
        X = np.array(X)
        Y = np.array(Y)
        peplen = np.array(peplen)
        print (X.shape, Y.shape)
        return X,Y,peplen
    
    # model
    def BuildModel(self, input_shape, nlayers = 1):
        print('BuildModel ... ')
        self.layer_count = nlayers
        self.model = Sequential()
        self.model.add(Masking(mask_value = 0, input_shape = input_shape))
        
        for i in range(nlayers):
            self.model.add(Bidirectional(LSTM(self.layer_size, return_sequences = True, dropout_W = 0.2, dropout_U = 0.2)))
            # self.model.add(LSTM(self.layer_size, return_sequences = True))
            self.model.add(Activation('relu'))
            self.model.add(Dropout(0.3))
            
            # self.model.add(TimeDistributed(Dense(self.layer_size)))
            # self.model.add(Activation('relu'))
            # self.model.add(Dropout(0.3))
        
        # self.model.add(Bidirectional(LSTM(self.layer_size, return_sequences = True, dropout_W = 0.2, dropout_U = 0.2)))
        # self.model.add(Activation('relu'))
        # self.model.add(Dropout(0.3))
        # self.layer_count += 1
        
        # self.model.add(TimeDistributed(Dense(self.layer_size)))
        # self.model.add(Activation('relu'))
        # self.model.add(Dropout(0.3))
        # self.layer_count += 1
        
        #the output layer
        self.model.add(TimeDistributed(Dense(self.max_charge * len(self.ion_types))))
        self.model.add(Activation('relu'))
        
    def InitModel(self, optimizer = 'adam', init_weights_file = None):
        if optimizer == 'sgd':
            optimizer = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        
        # in keras/objectives.py
        # mean_squared_error / mse
        # root_mean_squared_error / rmse
        # mean_absolute_error / mae
        self.model.compile(loss='mae', optimizer=optimizer)
        if not init_weights_file is None:
            self.model.load_weights(init_weights_file)
    
    def TrainModel(self, X_train, y_train, valid = None):
        self.model.fit(X_train, y_train, batch_size=self.batch_size, nb_epoch=self.epoches+self.initial_epoch, validation_split=0.0, initial_epoch=self.initial_epoch)

    def SaveWeights(self, weights_file):
        self.model.save_weights(weights_file)
    # model
    def TestModel(self, X_test, y_test, peplen_test, X_train = None, y_train = None, peplen_train = None):
        toplot = []
        if X_train is not None:
            print("Training set sim:")
            predictions = self.ModelPredict(X_train)
            compRNN = self.ComparePredReal_RNN(predictions, y_train, peplen_train)
            toplot.append(compRNN[0])
            ResTrain = compRNN[1]
        else:
            ResTrain = (0,0,0,0)
        print("Test set sim:")
        predictions = self.ModelPredict(X_test)
        compRNN = self.ComparePredReal_RNN(predictions, y_test, peplen_test)
        toplot.append(compRNN[0])
        ResTest = compRNN[1]
        plt.figure()
        plt.boxplot(toplot,sym='')
        plt.savefig("./models/lstm_layer=%d-epoch=%d-ch=%d-type=%s.png" %(self.layer_count, self.epoches+self.initial_epoch, self.max_charge, "".join(self.ion_types)))
        return (ResTrain,ResTest)

    # model
    def Predict(self, X):
        predictions = self.model.predict(X, batch_size = self.batch_size)
        predictions[predictions < 0] = 0
        predictions[predictions > 1] = 1
        return predictions

    # model
    def SaveModel(self, filepath = None):
        if filepath is None:
            self.model.save("./models/lstm_layer=%d-epoch=%d-ch=%d-type=%s.lstm" %(self.layer_count, self.epoches+self.initial_epoch, self.max_charge, "".join(self.ion_types)))
        else:
            self.model.save(filepath)
    
    # model
    def LoadModel(self, modelfile):
        self.model = load_model(modelfile)

    # data
    def LoadDataFolder(self, dataset_folder):
        print("Loading data from folder...")
        filenames = []
        for input_file in os.listdir(dataset_folder):
            if fnmatch.fnmatch(input_file, "*.plabel"):
                filenames.append(os.path.join(dataset_folder, input_file))
        return self.DirectIon2Vector(filenames)

    # data
    def LoadDataFiles(self, filenames):
        print("Loading data from files...")
        return self.DirectIon2Vector(filenames)
    
if __name__ == "__main__":
    IonLSTM().Run()

