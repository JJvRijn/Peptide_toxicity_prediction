"""
This is a boilerplate pipeline 'data_processing'
generated using Kedro 0.18.7

Created on Sat Apr 15 17:38:00 2023

@author: JJvRijn
"""
import pandas as pd
from modlamp.descriptors import GlobalDescriptor
import re
import numpy as np
def Scraping_DBAASP (input_excel_location,cell_types ="human", remove_unnatural=True) -> pd.DataFrame:
    #import the excel, this can be retrieved from https://dbaasp.org/search by pressing download the dataset
    DBAASP_data = input_excel_location
    #extract the usefull columns
    DBAASP_usefull = DBAASP_data[["ID", "SEQUENCE", "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE",
                                  "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION",
                                 "HEMOLITIC CYTOTOXIC ACTIVITY - UNIT", 'HEMOLITIC CYTOTOXIC ACTIVITY - TARGET CELL' 
                                  ]]
    #make the data workable
    #by removing unknown values and empty values
    DBAASP_usefull = DBAASP_usefull.dropna(thresh=3)
    DBAASP_usefull.reset_index(drop=True)
    #removing the unselected cell lines
    human = []
    for j in DBAASP_usefull.index:
        if cell_types in DBAASP_usefull[
        'HEMOLITIC CYTOTOXIC ACTIVITY - TARGET CELL'][j].lower():
            human.append(j)
    human = list( dict.fromkeys(human) )
    DBAASP_usefull = DBAASP_usefull.loc[human]
    DBAASP_usefull.reset_index(drop=True)
    #MEC isnt universally defined so is removed as correct input is near impossible
    #without processing the indvidual papers
    DBAASP_usefull = DBAASP_usefull[DBAASP_usefull[
        "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"] != 'MEC']
    DBAASP_usefull = DBAASP_usefull[DBAASP_usefull[
        "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"] != '-']
    DBAASP_usefull = DBAASP_usefull[DBAASP_usefull[
        "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"] != 'nan']
    DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"
                  ] = DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"].replace([
        "EC50", "IC50", 'LD50', 'CC50', 'LC50', 'MHC', 'LC90', 'IC90', 'ED50', 'LD90']
        , ["50%",'50%', '50%', '50%', '50%', '15%', '90%', '90%', '50%', '90%'])
    #removes rows with NA in specific cols
    DBAASP_usefull = DBAASP_usefull.dropna(subset=["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"])
    DBAASP_usefull = DBAASP_usefull.dropna(subset=["HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"])
    DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
    #Extracting int for Lysis_value
    mask3 = []
    for i in DBAASP_usefull.index:
        try:
            if "-" in str(DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"][i]):
                DBAASP_usefull.loc[i, "hem_activity"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                    "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"][i]))[1]       
            else:
                DBAASP_usefull.loc[i, "hem_activity"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                    "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"][i]))[0]
        except:
            mask3.append(i)
    DBAASP_usefull = DBAASP_usefull.drop(mask3)
    DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
    if remove_unnatural is True:
        #remove unnatural Amino acids
        NAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S",
               "T", "W", "Y", "V", ]
        unvalid =[]
        for j in DBAASP_usefull.index:
            if (not all ([NAA_seq in NAA for NAA_seq in DBAASP_usefull.SEQUENCE[j]])):
                unvalid.append(j)
        DBAASP_usefull = DBAASP_usefull.drop(unvalid)
        DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
        #add molweight to covert ug/ml to uM
        glob = GlobalDescriptor(DBAASP_usefull.SEQUENCE.tolist())
        glob.calculate_MW(amide=True)
        DBAASP_usefull['MW'] = glob.descriptor
        #see what concentrations are in µg/ml instead of µM
        mask = DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - UNIT"] == 'µg/ml'
        #extract numerical to remove >, <, <= etc.
        DBAASP_usefull["CONCENTRATION_µM"] = 0
        for o in range(len(DBAASP_usefull)):
            if "-" in str(DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"][o]):
                 DBAASP_usefull.loc[o, "conc_num"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"][o]))[1]
            else:
                DBAASP_usefull.loc[o, "conc_num"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                    "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"][o]))[0]
        #calculate the µM
        DBAASP_usefull["CONCENTRATION_µM"] = 0
        DBAASP_usefull.loc[mask, "CONCENTRATION_µM"] = ((DBAASP_usefull[
            "conc_num"][mask]).astype(float)*1000) / DBAASP_usefull[mask]['MW']
        #fill in values already in µM
        mask2 = ~mask
        DBAASP_usefull.loc[mask2, "CONCENTRATION_µM"] = DBAASP_usefull["conc_num"][mask2].astype(float)
        DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
        DBAASP_usefull = DBAASP_usefull.rename(columns ={"Lysis_value":"hem_activity", "SEQUENCE":"seq"})
        DBAASP_final = DBAASP_usefull[["ID", "seq", "hem_activity", "CONCENTRATION_µM"]]
        return DBAASP_final
    else:
        DBAASP_usefull = DBAASP_usefull.rename(columns = {"Lysis_value":"hem_activity", "SEQUENCE":"seq"})
        DBAASP_final = DBAASP_usefull[["ID", "seq", "hem_activity", "conc_num", "HEMOLITIC CYTOTOXIC ACTIVITY - UNIT"]]
        return DBAASP_final

def Scraping_DRAMP (input_excel_location, remove_unnatural=True,) -> pd.DataFrame:
    #import the data
    DRAMP_data = input_excel_location
    #extracting the usefull data
    seq = DRAMP_data.Sequence
    toxic = DRAMP_data.Hemolytic_activity
    ID = DRAMP_data.DRAMP_ID
    df_DRAMP = pd.DataFrame([ID, seq, toxic], index = ['ID', 'seq', 'toxic'])
    df_DRAMP = df_DRAMP.transpose()
    df_DRAMP.seq = list(map(lambda x: x.upper(),df_DRAMP.seq))
    #dropping unknown activities
    df_DRAMP = df_DRAMP[df_DRAMP["toxic"] != "Not found"]
    df_DRAMP = df_DRAMP[df_DRAMP["toxic"] != " No hemolytic activity information found."]
    df_DRAMP = df_DRAMP[df_DRAMP["toxic"] != ' Not found in the literature']
    df_DRAMP = df_DRAMP.dropna(subset=['toxic'])
    df_DRAMP.toxic=df_DRAMP.toxic.str.replace(r'\[[^]]*\]', '')
    df_DRAMP = df_DRAMP.reset_index(drop=True)
    #extracting activities
    df_DRAMP['toxic2'] = 'NaN'
    df_DRAMP['toxic3'] = 'NaN'
    df_DRAMP['toxic4'] = 'NaN'
    new_values = pd.DataFrame(columns = df_DRAMP.columns)
    for n in df_DRAMP.index:
        df_DRAMP.toxic[n] = re.sub('\[.*?\]', '', df_DRAMP.toxic[n])
        if "-" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic']= re.sub('-[0-9. ]+[%%μµnm ]', '', df_DRAMP.toxic[n])
        if "±" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic']= re.sub('±[0-9. ]+[%%μµnm ]', '', df_DRAMP.toxic[n])    
        if "no hemolytic activity" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('no hemolytic activity', '0% hem')
        if "no detectable hemo" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('no detectable '
                                                                 , '0%')
        if "MHC" in df_DRAMP.toxic[n]:
            if "MHC10" not in df_DRAMP.toxic[n]:
                df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('MHC', '15% hem')
        if "Non-hemoly" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('Non-', '0%')
        
        df_DRAMP.toxic2[n] = re.findall(r"[-+]?\d*\.\d+|\d+", str(df_DRAMP.toxic[n]))
        if (len(df_DRAMP.toxic2[n])==2):
            correct_row = df_DRAMP.loc[n]
            correct_row['toxic3'] = (df_DRAMP.toxic2[n])[0]
            correct_row['toxic4'] = (df_DRAMP.toxic2[n])[1]
            new_values.loc[len(new_values)] = correct_row            
        if (len(df_DRAMP.toxic2[n])>2):
            if (len(df_DRAMP.toxic2[n])%2 == 0):
                single_pep = df_DRAMP.toxic[n].replace('and',',').replace('##',',').replace(';',',').split(sep =",")
                for r in single_pep:
                    two_things = re.findall(r"[-+]?\d*\.\d+|\d+", str(r))                
                    if len(two_things) == 2:
                        if "positive control" not in r:
                            correct_row = df_DRAMP.loc[n]
                            correct_row['toxic'] = r
                            correct_row['toxic2'] = two_things
                            correct_row['toxic3'] = two_things[0]
                            correct_row['toxic4'] = two_things[1]
                            new_values.loc[len(new_values)] = correct_row
                    if len(two_things) == 4:
                        if two_things[2] == "10":
                            correct_row = df_DRAMP.loc[n]
                            correct_row['toxic'] = r
                            correct_row['toxic2'] = two_things
                            correct_row['toxic3'] = (df_DRAMP.toxic2[n])[0]
                            correct_row['toxic4'] = float((df_DRAMP.toxic2[n])[1])*10**float((df_DRAMP.toxic2[n])[3])
                            new_values.loc[len(new_values)] = correct_row
    df_DRAMP = new_values
    #extract the unit of concentration
    df_DRAMP['toxic5'] = 'NaN'
    for n in range(len(df_DRAMP)):
        if "M" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "M"
        if "g/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "g/ml"
        if "g/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mg/ml" 
        if "mol/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "M"   
        if "μM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"
        if "µM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"
        if "µmol/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"
        if "μmol/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"        
        if "µg/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"
        if "μg/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"
        if "μg/mL" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"
        if "µg/mL" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"        
        if "nM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "nM"       
        if "mg/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mg/L"
        if "mg/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mg/ml"
        if "ng/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "ng/ml" 
        if "nmol/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM" 
        if "mM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mM"              
    df_DRAMP= df_DRAMP[df_DRAMP.toxic5 != "NaN"]
    df_DRAMP = df_DRAMP.reset_index(drop=True)
    #remove the unnatural AA or return current info
    if remove_unnatural is True:
        NAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S",
               "T", "W", "Y", "V", ]
        unvalid =[]
        for j in df_DRAMP.index:
            if (not all ([NAA_seq in NAA for NAA_seq in df_DRAMP.seq[j]])):
                unvalid.append(j)
        df_DRAMP = df_DRAMP.drop(unvalid)
        df_DRAMP = df_DRAMP.reset_index(drop=True)
    else:
         df_final = df_DRAMP[["ID", "seq", "toxic3", "toxic4", "toxic5"]]
         df_final = df_final.rename( columns={"toxic3":"hem_activity", "toxic4":"concentration", "toxic5":"unit_concentration"})
         return df_final
    glob = GlobalDescriptor((df_DRAMP.seq.tolist()))
    glob.calculate_MW(amide=True)
    df_DRAMP['MW'] = glob.descriptor
    #see what concentrations are not in µM
    mask0 = df_DRAMP["toxic5"] == 'μg/ml'
    mask1 = df_DRAMP["toxic5"] == 'mg/L'
    mask2 = df_DRAMP["toxic5"] == 'mg/ml'
    mask3 = df_DRAMP["toxic5"] == 'μM'
    mask4 = df_DRAMP["toxic5"] == 'nM'
    mask5 = df_DRAMP["toxic5"] == 'mM'
    mask6 = df_DRAMP["toxic5"] == 'M'
    mask7 = df_DRAMP["toxic5"] == 'g/ml'
    mask8 = df_DRAMP["toxic5"] == 'ng/ml'
    #calculate the µM
    df_DRAMP["CONCENTRATION_µM"] = "NaN"
    df_DRAMP.loc[mask0, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask0]).astype(float)*1000) / df_DRAMP[mask0]['MW']
    df_DRAMP.loc[mask1, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask1]).astype(float)*1000) / df_DRAMP[mask1]['MW']
    df_DRAMP.loc[mask2, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask2]).astype(float)*1000000) / df_DRAMP[mask2]['MW']
    df_DRAMP.loc[mask4, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask4]).astype(float)/1000)
    df_DRAMP.loc[mask5, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask5]).astype(float)*1000)
    df_DRAMP.loc[mask6, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask6]).astype(float)*1000000)
    df_DRAMP.loc[mask7, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask7]).astype(float)*1000000000) / df_DRAMP[mask7]['MW']
    df_DRAMP.loc[mask8, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask8]).astype(float)) / df_DRAMP[mask8]['MW']
    #fill in values already in µM
    df_DRAMP.loc[mask3, "CONCENTRATION_µM"] = df_DRAMP["toxic4"][mask3].astype(float)
    df_final = df_DRAMP[["ID", "seq", "toxic3", "CONCENTRATION_µM"]]
    df_final = df_final.rename(columns={"toxic3":"hem_activity"})
    return df_final


def string_to_act_conc (df_hemolytik) -> pd.DataFrame:
    #input is a dataframe with a column that is called toxic containing the string that explains the toxicity    
    #automatic reading of the activity from string to floats
    df_hemolytik =  df_hemolytik.rename(columns={"hemolycity string" : "toxic", "SEQUENCE":"seq"})
    df_hemolytik['toxic2'] = np.nan
    df_hemolytik['toxic3'] = "other"
    df_hemolytik['toxic4'] = "other"
    df_hemolytik['toxic5'] = "other"
    for n in range(len(df_hemolytik)):
        if "-" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic']= re.sub('-[^>]+[%μµnm ]', '', df_hemolytik.toxic[n])
        if "±" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic']= re.sub('±[^>]+[%μµnm ]', '', df_hemolytik.toxic[n])
        if "MHC" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('MHC', '15% hem')
        if "No hemolysis" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('No hemolysis', '0% hem')
        if "Lethal" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('Lethal', '100% hem')
        if "not detected" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('not detected', '0% hem')
        if "M" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "M"
        if "μM" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "μM"
        if "µM" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "μM"
        if "µg/ml" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "μg/ml"
        if "μg/ml" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "μg/ml"
        if "μg/mL" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "μg/ml"
        if "µg/mL" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "μg/ml"        
        if "nM" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "nM"
        if "mg/L" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "mg/L"
        if "mg/ml" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "mg/ml"        
        if "mM" in df_hemolytik.toxic[n]:
            df_hemolytik.loc[n, 'toxic3'] = "mM"
        df_hemolytik.toxic2[n] = re.findall(r"[-+]?\d*\.\d+|\d+", str(df_hemolytik.toxic[n]))
        if (len(df_hemolytik.toxic2[n])==2):
            df_hemolytik.loc[n, 'toxic4'] = (df_hemolytik.toxic2[n])[0]
            df_hemolytik.loc[n, 'toxic5'] = (df_hemolytik.toxic2[n])[1]
    #removing peptides that could not have their activity added correctly
    unvalid1= []
    unvalid2= []
    for n in range(len(df_hemolytik)):
        if ((df_hemolytik.toxic3[n]== 'other')):
            unvalid1.append(n)
        if ((df_hemolytik.toxic4[n]== 'other')):
            unvalid2.append(n)
    unvalid = unvalid1 + list(set(unvalid2) - set(unvalid1))
    df_hemolytik = df_hemolytik.drop(unvalid)
    df_hemolytik.reset_index(drop=True)        
    #standardize the activity
    df_hemolytik = df_hemolytik.reset_index(drop=True)
    #types of units that need conversion "μg/ml", "nM", "M", "mM", "mg/L", "mg/ml"
    glob = GlobalDescriptor(df_hemolytik.seq.tolist())
    glob.calculate_MW(amide=True)
    df_hemolytik['MW'] = glob.descriptor
    #see what concentrations are in µg/ml instead of µM
    mask0 = df_hemolytik["toxic3"] == 'µg/ml'
    mask1 = df_hemolytik["toxic3"] == 'mg/L'
    mask2 = df_hemolytik["toxic3"] == 'mg/ml'
    mask3 = df_hemolytik["toxic3"] == 'μM'
    mask4 = df_hemolytik["toxic3"] == 'nM'
    mask5 = df_hemolytik["toxic3"] == 'mM'
    mask6 = df_hemolytik["toxic3"] == 'M'
    #calculate the µM
    df_hemolytik["CONCENTRATION_µM"] = 0
    df_hemolytik.loc[mask0, "CONCENTRATION_µM"] = (
        (df_hemolytik["toxic5"][mask0])*1000) / df_hemolytik[mask0]['MW']
    df_hemolytik.loc[mask1, "CONCENTRATION_µM"] = (
        (df_hemolytik["toxic5"][mask1]).astype(float)*1000) / df_hemolytik[mask1]['MW']
    df_hemolytik.loc[mask2, "CONCENTRATION_µM"] = (
        (df_hemolytik["toxic5"][mask2]).astype(float)*1000000) / df_hemolytik[mask2]['MW']
    df_hemolytik.loc[mask4, "CONCENTRATION_µM"] = (
        (df_hemolytik["toxic5"][mask4]).astype(float)/1000)
    df_hemolytik.loc[mask5, "CONCENTRATION_µM"] = (
        (df_hemolytik["toxic5"][mask5]).astype(float)*1000)
    df_hemolytik.loc[mask6, "CONCENTRATION_µM"] = (
        (df_hemolytik["toxic5"][mask6]).astype(float)*1000000)
    #fill in values already in µM
    df_hemolytik.loc[mask3, "CONCENTRATION_µM"] = df_hemolytik["toxic5"][mask3].astype(float)
    df_hemolytik.drop(columns = ["toxic2", "toxic3", "toxic5", "MW"])
    df_hemolytik.rename(columns={"toxic4": "hem_activity"})
    return df_hemolytik

def get_var_name(variable):
     for name, value in globals().items():
        if value is variable:
            return name

def concat_w_source(*args )-> pd.DataFrame:
    import pandas as pd
    objs = list(args)
    b = [get_var_name(el) for el in objs]
    for i in range(len(objs)):
        objs[i]["source"] = b[i]
    concatted = pd.concat(objs, ignore_index=True).astype({'ID':str, 'CONCENTRATION_µM':float, "hem_activity":float})
    return concatted[["ID", "seq", "hem_activity", 	"CONCENTRATION_µM",]]
