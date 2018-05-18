from .AAMass import AAMass

aamass = AAMass()

def AAMass_fix_C57():
    aamass.aa_mass_dict['C'] += 57.021464

def calc_mod_mass_list(peptide, modinfo):
    items = modinfo.split(";")
    modlist = []
    for mod in items:
        if mod != '':
            site, modname = mod.split(",")
            site = int(site)
            modlist.append( (site, modname) )
    modlist.sort()
    
    modmass = [0]*(len(peptide)+2)
    lossmass = [0]*(len(peptide)+2)
    for mod in modlist:
        modmass[mod[0]] = aamass.mod_mass_dict[mod[1]][0]
        lossmass[mod[0]] = aamass.mod_mass_dict[mod[1]][1]
    return modmass,lossmass
    
def calc_b_ions(peptide, modinfo):
    modmass_list, modloss_list = calc_mod_mass_list(peptide, modinfo)
    b_ions = []
    mass_nterm = modmass_list[0]
    for i in range(len(peptide)-1):
        mass_nterm += aamass.aa_mass_dict[peptide[i]] + modmass_list[i+1]
        b_ions.append(mass_nterm)
    pepmass = b_ions[-1] + aamass.aa_mass_dict[peptide[-1]] + modmass_list[len(peptide)] + modmass_list[len(peptide)+1] + aamass.mass_H2O
    return b_ions, pepmass
    
def calc_y_from_b(bions, pepmass):
    return [pepmass - b for b in bions]
    
def calc_c_from_b(bions, pepmass = 0):
    return [b + aamass.mass_NH3 for b in bions]
    
def calc_z_from_b(bions, pepmass):
    return [pepmass - b - aamass.mass_NH3 + aamass.mass_H for b in bions]
    