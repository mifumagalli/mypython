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

        dum=fle.readline()
        hd_name=fle.readline()
        dum=fle.readline()
        hd_isf=fle.readline()
        dum=fle.readline()
        hd_toma=fle.readline()
        dum=fle.readline()
        hd_sfr=fle.readline()
        dum=fle.readline()
        hd_ninterv=fle.readline()
        dum=fle.readline()
        hd_xponent=fle.readline()
        dum=fle.readline()
        hd_xmaslim=fle.readline()
        dum=fle.readline()
        hd_sncut=fle.readline()
        dum=fle.readline()
        hd_bhcut=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        hd_iz=fle.readline()
        dum=fle.readline()
        hd_iwind=fle.readline()
        dum=fle.readline()
        hd_time1=fle.readline()
        dum=fle.readline()
        hd_jtime=fle.readline()
        dum=fle.readline()
        hd_tbiv=fle.readline()
        dum=fle.readline()
        hd_itbiv=fle.readline()
        dum=fle.readline()
        hd_tmax=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        hd_jmg=fle.readline()
        dum=fle.readline()
        hd_lminmax=fle.readline()
        dum=fle.readline()
        hd_tdel=fle.readline()
        dum=fle.readline()
        hd_iatmos=fle.readline()
        dum=fle.readline()
        dum=fle.readline()
        hd_ilib=fle.readline()
        dum=fle.readline()
        hd_iline=fle.readline()
        dum=fle.readline()
        hd_ivtrsg=fle.readline()
        dum=fle.readline()
        hd_io=fle.readline()
        fle.close()

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
        if(files[0] == '1'):
            
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
        if(files[6] == '1'):
            
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

                
