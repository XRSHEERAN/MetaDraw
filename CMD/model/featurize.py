import math
import numpy as np

from .fragmentation_config import Common_Config

class Ion2Vector:
    def __init__(self, modini=Common_Config().modini, ion_types = ["y"], max_charge = 1, cleave_site_num = 19, prev = 1, next = 1, from_CTerm = False, fixmod = ["Carbamidomethyl[C]"], varmod = []):
        self.ion_types = ion_types
        self.ion_terms = Common_Config().ion_terms
        self.max_charge = max_charge
        self.ionlen = cleave_site_num # ionlen is the number of cleavage sites
        self.prev = prev
        self.next = next
        self.from_CTerm = from_CTerm
        self.fixmod = fixmod
        self.varmod = varmod
        self.aa2vector = self.AAVectorDict()
        self.AA_idx = dict(zip("ACDEFGHIKLMNPQRSTVWY",range(0,len(self.aa2vector))))
        self.modini = modini
        self.__read_mod__()
        self.__get_ion_type_init__()
         
    def __read_mod__(self):
        # feature_vector of mod: [#C, #H, #N, #O, #S, #P, #metal, #other]
        common = ['C','H','N','O','S','P']
        metal = ['Na','Ca','Fe','K','Mg']
        self.mod_feature_size = len(common) + 2
        f = open(self.modini)
        lines = f.readlines()
        f.close()
        
        def parse_element(elem_str):
            feature = [0] * self.mod_feature_size
            elems = elem_str.split(')')[:-1]
            for elem in elems:
                chem, num = elem.split('(')
                num = int(num)
                
                try:
                    idx = common.index(chem)
                    feature[idx] = num
                except:
                    if chem in metal:
                        feature[-2] += num
                    else:
                        feature[-1] += num
            return np.array(feature, dtype=np.float32)
        
        self.mod_feature = {}
        num_mod = int(lines[0].split("=")[-1])
        for i in range(num_mod):
            modinfo = lines[i*2 + 2].split(' ')
            modname = modinfo[0].split('=')[0]
            elem_str = modinfo[-1].strip()
            self.mod_feature[modname] = parse_element(elem_str)
    
    def __get_ion_type_init__(self):
        self.sorted_ion_types = sorted(self.ion_types, key = lambda x: -len(x))
    def __get_ion_type__(self, ion_type):
        for type in self.sorted_ion_types:
            if type in ion_type:
                return type
        return ""

    def AAVectorDict(self):
        aa2vector_map = {}
        s = "ACDEFGHIKLMNPQRSTVWY"
        v = [0]*len(s)
        v[0] = 1
        for i in range(len(s)):
            aa2vector_map[s[i]] = list(v)
            v[i],v[(i+1) % 20] = 0,1
        return aa2vector_map

    def Parse_Ion2vector(self, peptide, ionidx, charge):
        seqidx = ionidx
        # look at the ion's previous "prev" N_term AA
        v = []
        for i in range(seqidx - self.prev, seqidx):
            if i < 0:
                v.extend([0]*len(self.aa2vector))
            else:
                v.extend(self.aa2vector[peptide[i]])
        # look at the ion's next "next" C_term AAs
        for i in range(seqidx, seqidx + self.next):
            if i >= len(peptide):
                v.extend([0]*len(self.aa2vector))
            else:
                v.extend(self.aa2vector[peptide[i]])
        
        #the number of each AA before "prev" in NTerm
        NTerm_AA_Count = [0]*len(self.aa2vector)
        for i in range(seqidx - self.prev):
            NTerm_AA_Count[self.AA_idx[peptide[i]]] += 1
        v.extend(NTerm_AA_Count)
        
        #the number of each AA after "next" in CTerm
        CTerm_AA_Count = [0]*len(self.aa2vector)
        for i in range(seqidx + self.next, len(peptide)):
            CTerm_AA_Count[self.AA_idx[peptide[i]]] += 1
        v.extend(CTerm_AA_Count)
        
        if ionidx == 1: CTerm = 1
        else: CTerm = 0
        if ionidx == len(peptide)-1: NTerm = 1
        else: NTerm = 0
        v.extend([NTerm,CTerm])
        
        vchg = [0]*6
        vchg[charge-1] = 1
        v.extend(vchg)
        return np.array(v, dtype=np.float32)
        
    def PaddingZero(self, ionvec, intenvec, modvec):
        for i in range(self.ionlen - len(ionvec)):
            ionvec.append(np.array([0] * ( (len(ionvec[0]) ) ), dtype=np.float32))
            intenvec.append(np.array([0] * (self.max_charge * len(self.ion_types)), dtype=np.float32))
            modvec.append(np.array([0] * ( (len(modvec[0]) ) ), dtype=np.float32))
        return (ionvec, intenvec, modvec)

    def FeaturizeOnePeptide(self, peptide, modinfo, charge):
        if len(peptide) > self.ionlen + 1: return None, None
        if charge > 6: return None, None
        modlist = modinfo.split(';')
        mod_idx_feature = {}
        for mod in modlist:
            if not mod: continue
            modtmp = mod.split(',')
            idx = int(modtmp[0])
            modname = modtmp[1]
            if modname in self.fixmod:
                continue
            else:
                if idx ==0 or idx == 1: idx == 0
                elif idx >= len(peptide): idx = len(peptide)-1
                else: idx -= 1
                mod_idx_feature[idx] = self.mod_feature[modname]
                
        x = []
        mod_x = []
        for ionidx in range(1, len(peptide)):
            if ionidx-1 in mod_idx_feature:
                mod_x.append(mod_idx_feature[ionidx-1])
            else:
                mod_x.append(np.array([0]*self.mod_feature_size, dtype=np.float32))
                    
            v = self.Parse_Ion2vector(peptide, ionidx, charge)
            x.append(v)
        x, __, mod_x = self.PaddingZero(x, [], mod_x)
        # for i in range(len(ionvec)):
            # out1.write("%f,%s\n" %(intenvec[i], ionvec[i]))
        if self.from_CTerm:
            x = x[::-1] # from C-Term to N-Term
            mod_x = mod_x[::-1]
        return x, mod_x
    
    def Featurize_ignore_mod(self, ion_file, RNN = True):
        f = open(ion_file)
        
        X = []
        Y = []
        peplen = []

        while True:
            line = f.readline()
            if line == "": break
            type2inten = {}
            allinten = []
            items = line.split("\t")
            peptide = items[1]
            if len(peptide) > self.ionlen + 1: continue
            charge = int(items[0].split(".")[-3])
            if charge > 6: continue
            ion_types = items[-2].split(",")[:-1]
            if len(ion_types) < len(peptide): continue
            
            intens = [float(item) for item in items[-1].split(",")[:-1]]
            for i in range(len(ion_types)):
                ion_type = self.__get_ion_type__(ion_types[i])
                if ion_type in self.ion_types:
                    allinten.append(intens[i])
                if not ion_type in self.ion_types: continue
                if int(ion_types[i][ion_types[i].find("+")+1:]) > self.max_charge: continue
                type2inten[ion_types[i]] = intens[i]
            
            if len(type2inten) < len(peptide) / 2: continue
            
            maxinten = max(allinten)
            
            ionvec = []
            intenvec = []
                    
            for ionidx in range(1, len(peptide)):
                        
                v = self.Parse_Ion2vector(peptide, ionidx, charge)
                ionvec.append(v)
                v = []
                for ion_name in self.ion_types:
                    for ch in range(1, self.max_charge + 1):
                        if self.ion_terms[ion_name] == 'c':
                            iontype = "{}{}+{}".format(ion_name, len(peptide) - ionidx, ch)
                        else:
                            iontype = "{}{}+{}".format(ion_name, ionidx, ch)
                        if iontype in type2inten:
                            v.append(type2inten[iontype] / maxinten)
                        else:
                            v.append(0)
                intenvec.append(np.array(v, dtype=np.float32))
            
            # for i in range(len(ionvec)):
                # out1.write("%f,%s\n" %(intenvec[i], ionvec[i]))
            if self.from_CTerm:
                ionvec = ionvec[::-1] # from C-Term to N-Term
                intenvec = intenvec[::-1] # from C-Term to N-Term
            #################### RNN Feature ###############################
            if RNN:
                ionvec, intenvec, __ = self.PaddingZero(ionvec, intenvec,[])
                X.append(ionvec)
                Y.append(intenvec)
            else:
                X.extend(ionvec)
                Y.extend(intenvec)
            peplen.append(len(peptide))

        f.close()
        return (X,Y,peplen)
    
    def Featurize(self, ion_file, RNN = True):
        f = open(ion_file)
        
        X = []
        Y = []
        peplen = []
        
        if len(self.varmod) > 0: mod_X = []

        while True:
            line = f.readline()
            if line == "": break
            type2inten = {}
            allinten = []
            items = line.split("\t")
            peptide = items[1]
            if len(peptide) > self.ionlen + 1: continue
            charge = int(items[0].split(".")[-3])
            if charge > 6: continue
            ion_types = items[-2].split(",")[:-1]
            if len(ion_types) < len(peptide): continue
            # print(items[0])
            modlist = items[2].split(';')
            mod_idx_feature = {}
            if len(self.varmod) > 0 and len(modlist) == 0: continue
            elif len(modlist) != 0 and len(self.varmod) > 0:
                unexpected_mod = False
                for mod in modlist:
                    modtmp = mod.split(',')
                    idx = int(modtmp[0])
                    modname = modtmp[1]
                    if modname in self.fixmod:
                        continue
                    elif modname in self.varmod:
                        if idx ==0 or idx == 1: idx == 0
                        elif idx >= len(peptide): idx = len(peptide)-1
                        else: idx -= 1
                        mod_idx_feature[idx] = self.mod_feature[modname]
                    else:
                        unexpected_mod = True
                        break
                if unexpected_mod: continue
            
            intens = [float(item) for item in items[-1].split(",")[:-1]]
            for i in range(len(ion_types)):
                ion_type = self.__get_ion_type__(ion_types[i])
                if ion_type in self.ion_types:
                    allinten.append(intens[i])
                if not ion_type in self.ion_types: continue
                if int(ion_types[i][ion_types[i].find("+")+1:]) > self.max_charge: continue
                type2inten[ion_types[i]] = intens[i]
            
            if len(type2inten) < len(peptide) / 2: continue
            
            maxinten = max(allinten)
            
            ionvec = []
            intenvec = []
            modvec = []
                    
            for ionidx in range(1, len(peptide)):
                if ionidx-1 in mod_idx_feature:
                    modvec.append(mod_idx_feature[ionidx-1])
                else:
                    modvec.append(np.array([0]*self.mod_feature_size, dtype=np.float32))
                        
                v = self.Parse_Ion2vector(peptide, ionidx, charge)
                ionvec.append(v)
                v = []
                for ion_name in self.ion_types:
                    for ch in range(1, self.max_charge + 1):
                        if self.ion_terms[ion_name] == 'c':
                            iontype = "{}{}+{}".format(ion_name, len(peptide) - ionidx, ch)
                        else:
                            iontype = "{}{}+{}".format(ion_name, ionidx, ch)
                        if iontype in type2inten:
                            v.append(type2inten[iontype] / maxinten)
                        else:
                            v.append(0)
                intenvec.append(np.array(v, dtype=np.float32))
            
            # for i in range(len(ionvec)):
                # out1.write("%f,%s\n" %(intenvec[i], ionvec[i]))
            if self.from_CTerm:
                ionvec = ionvec[::-1] # from C-Term to N-Term
                intenvec = intenvec[::-1] # from C-Term to N-Term
                modvec = modvec[::-1]
            #################### RNN Feature ###############################
            if RNN:
                ionvec, intenvec, modvec = self.PaddingZero(ionvec, intenvec,modvec)
                X.append(ionvec)
                Y.append(intenvec)
                if len(self.varmod) > 0: mod_X.append(modvec)
            else:
                X.extend(ionvec)
                Y.extend(intenvec)
                mod_X.extend(modvec)
            peplen.append(len(peptide))

        f.close()
        if len(self.varmod) == 0: return (X,Y,peplen)
        else: return (X,mod_X,Y,peplen)
        
if __name__ == '__main__':
    i2v = Ion2Vector()
    i2v.Featurize(r'f:\zhouxiexuan\Mann-Synthetic-NBT-2013-Velos-40-pf2\Phospho\1_HCDFT.pf2.HCD.plabel')
