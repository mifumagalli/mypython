class Sb99:
    
    "Define Class"
    def __init__(self):
        self.model=""
       
    def loadmod(self,name,path='./'):
        
        from astropy.table import Table
        import numpy as np

        #store model name
        self.model=name

        #parse the input file 
        fle=open(path+name+'.input1','r')
        alllines=[]
        for ll in fle:
            alllines.append(ll)
        fle.close

        for ii,ll in enumerate(alllines):
            if('[NAME]' in ll):
                hd_name=alllines[ii+1].rstrip()
            if('[ISF]' in ll):
                hd_isf=alllines[ii+1].rstrip()
            if('[TOMA]' in ll):
                hd_toma=alllines[ii+1].rstrip()
            if('[SFR]' in ll):
                hd_sfr=alllines[ii+1].rstrip()
            if('[NINTERV]' in ll):
                hd_ninterv=alllines[ii+1].rstrip()
            if('[XPONENT]' in ll):
                hd_xponent=alllines[ii+1].rstrip()
            if('[XMASLIM]' in ll):
                hd_xmaslim=alllines[ii+1].rstrip()
            if('[SNCUT]' in ll):
                hd_sncut=alllines[ii+1].rstrip()
            if('[BHCUT]' in ll):
                hd_bhcut=alllines[ii+1].rstrip()
            if('[IZ]' in ll):
                hd_iz=alllines[ii+7].rstrip()
            if('[IWIND]' in ll):
                hd_iwind=alllines[ii+1].rstrip()
            if('[TIME1]' in ll):
                hd_time1=alllines[ii+1].rstrip()
            if('[JTIME]' in ll):
                hd_jtime=alllines[ii+1].rstrip()
            if('[TBIV]' in ll):
                hd_tbiv=alllines[ii+1].rstrip()
            if('[ITBIV]' in ll):
                hd_itbiv=alllines[ii+1].rstrip()
            if('[TMAX]' in ll):
                hd_tmax=alllines[ii+1].rstrip()
            if('[JMG]' in ll):
                hd_jmg=alllines[ii+1].rstrip()
            if('[LMIN,' in ll):
                hd_lminmax=alllines[ii+1].rstrip()
            if('[TDEL]' in ll):
                hd_tdel=alllines[ii+1].rstrip()
            if('[IATMOS]' in ll):
                hd_iatmos=alllines[ii+1].rstrip()
            if('[ILIB]' in ll):
                hd_ilib=alllines[ii+2].rstrip()
            if('[ILINE]' in ll):
                hd_iline=alllines[ii+1].rstrip()
            if('[IVT,' in ll):
                hd_ivtrsg=alllines[ii+1].rstrip()
            if('[IO1,' in ll):
                hd_io=alllines[ii+1].rstrip()
        
        #store header
        self.header={"NAME":hd_name,"ISF":hd_isf,"TOMA":hd_toma,
                     "SFR":hd_sfr,"NINTERV":hd_ninterv,"XPONENT":hd_xponent,
                     "XMASLIM":hd_xmaslim,"SNCUT":hd_sncut,"BHCUT":hd_bhcut,
                     "IZ":hd_iz,"IWIND":hd_iwind,"TIME1":hd_time1,
                     "JTIME":hd_jtime,"TBIV":hd_tbiv,"ITBIV":hd_itbiv,
                     "TMAX":hd_tmax,"JMG":hd_jmg,"LMINMAX":hd_lminmax,
                     "TDEL":hd_tdel,"IATMOS":hd_iatmos,"ILIB":hd_ilib,
                     "ILINE":hd_iline,"IVTRSG":hd_ivtrsg,"IO":hd_io}
        
        #parse other files 
        files=(self.header["IO"].rstrip()).split(",")
        
        #quanta file
        if(files[0] == '+1'):
            
            fle=open(path+name+'.quanta1','r')
            
            time=[]
            HI=[]
            HeI=[]
            HeII=[]
            LogL=[]

            for line in fle:
                line = line.rstrip()
                columns = line.split()
                if(len(columns) == 8): 
                    time.append(float(columns[0]))
                    HI.append(float(columns[1]))
                    HeI.append(float(columns[3]))
                    HeII.append(float(columns[5]))
                    LogL.append(float(columns[7]))
            fle.close()
            self.quanta={"TIME":time,"HI":HI,"HeI":HeI,"HeII":HeII,"LogL":LogL}
        
        #spectra file
        if(files[6] == '+1'):
            
            fle=open(path+name+'.spectrum1','r')
            
            time=[]
            wavelength=[]
            total=[]
            stellar=[]
            nebular=[]

            for line in fle:
                line = line.rstrip()
                columns = line.split()
                if(len(columns) == 5): 
                    time.append(float(columns[0]))
                    wavelength.append(float(columns[1]))
                    total.append(float(columns[2]))
                    stellar.append(float(columns[3]))
                    nebular.append(float(columns[4]))
            fle.close()
            
            #find times and store in tabular form
            untimes=sorted(set(time)) 
            ntimes=len(untimes)
            
            #go to numpy arrays
            time=np.array(time)
            wavelength=np.array(wavelength)
            total=np.array(total)
            stellar=np.array(stellar)
            nebular=np.array(nebular)

            self.spectrum=[]
            for tt in untimes:
                nzero=np.nonzero( abs(time-tt) < 1)
                nzero=nzero[0]
                self.spectrum.append({"TIME":time[nzero],"WAVE":wavelength[nzero],"TOTAL":total[nzero],
                                      "STELLAR":stellar[nzero],"NEBULAR":nebular[nzero]})

                
