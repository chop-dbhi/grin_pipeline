class gender:
    def __init__(self,sampletable):
        self.sampletable = sampletable

    def getGender(self,proband):
        res = self.sampletable.loc[self.sampletable['Subject'] ==  proband,"Sex"]
        sex = str(res.values[0])
        return sex

    # trio - EG0013_EG0013P
    # sampletable - pandas table
    def trioGender(self,trio):
        proband=trio.split("_")[1]
        return self.getGender(proband)

    def getFamilyFromMember(self,member):
        rows = self.sampletable.loc[((self.sampletable['Subject'] == member) | (self.sampletable['Mother'] == member) | (self.sampletable['Father'] == member)) & ((self.sampletable['Mother'] != '') & (self.sampletable['Father'] != ''))]
        family = "{0}_{1}".format(rows['FamilyID'].tolist()[0],rows['Subject'].tolist()[0])
        return family

    #trio give as string_proband
    def getMembersFromTrio(self,trio):
        family,proband = trio.split("_")
        rows = self.sampletable.loc[(self.sampletable['FamilyID'] == family) & (self.sampletable['Subject'] == proband)]
        mom = rows['Mother'].tolist()[0]
        dad = rows['Father'].tolist()[0]
        pro = rows['Subject'].tolist()[0]
        return (pro,mom,dad)