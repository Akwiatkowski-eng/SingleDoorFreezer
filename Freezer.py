# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 11:15:06 2020

@author: Adrian
"""

import CoolProp.CoolProp as CP
import basicTermFlowLib as term
#import cabinet as cabinet
#import compressor as compressor
#import capillary as capillary
#import HEcalc as HE
import math as m

'''
Initial and boundary conditions
'''




T_env = 273+24                              #Temperature in room
T_inside = 273+20                        #Initial temperature inside appliance
#p_boil = 0.59*10**5                           #Initial low side pressure
T_inlet = 273                                #initial temperature at compressor inlet
mass_air = 0.05                              #mass flow through evaporator kg/s constant value for now
p_airin = 101300                            #pressure of air inside cabinet
p_crit = CP.PropsSI('PCRIT','R600a')        #critical refrigerant pressure - needed for evaporator
p_cond = 6*10**5                            #Initial condensation pressure
delta_time = 1                              #initial time step
time = 0                                    #start time of simulation
end_time = 8*3600                             #end time of simulation                            
T_stat = 273-16                             #desired compartment temperature
T_var1 = T_stat - 1                         #Actual temperature which stops the compressor
T_var2 = T_stat + 1                         #Actual temperature which starts the compressor
delta_t = 300

class fin_evap:
    
    HEtype = 'fined evaporator'
    p_crit = CP.PropsSI('PCRIT','R600a')
    p_r = 0.59*10**5/p_crit
    M = CP.PropsSI('MOLARMASS','R600A')*1000 #kg/mol
    
    def __init__(self,d_ex,d_in,pitch,tube_length,n_rows,n_columns,s_1,s_2, ):
        self.d_ex = d_ex*10**(-3)   
        self.d_in = d_in*10**(-3)                  
        self.pitch = pitch*10**(-3)                  
        
        self.tube_length = tube_length*10**(-3)
        self.n_rows = n_rows
        self.n_columns = n_columns
        self.s_1 = s_1*10**(-3)
        self.s_2 = s_2*10**(-3)
        self.fin_height = self.s_1/2
        self.Area_air = self.s_1*self.n_columns*self.tube_length       #cross section for air flow at inlet
        self.hydDiam_air = 2*(s_1*self.fin_height)**0.5      #hydraulic diameter for air
        self.Area_ref = 3.14*self.d_in**2/4
        self.lambda_wall = 220
        self.delta = 0.00015
        self.Area_ext = (2*(self.s_1*self.s_2-3.14/4*self.d_ex**2)+3.14*self.d_ex*(self.pitch-self.delta))*(self.tube_length/self.pitch)*self.n_rows*self.n_columns
        self.Area_in = 3.14*self.d_in*self.tube_length*self.n_columns*self.n_rows
        
    #def rho(self, s_2, d_ex):
        #rho = s_2/d_ex
        #return rho
        
#evap = fin_evap(5,6,10,400,10,2,10,15)
#tho = evap.rho(evap.s_2,evap.d_ex)
#print(evap.d_ex,evap.s_2,tho)
        
    def Nusselt_air(self,d_ex,pitch,fin_height,Re):
        c = 0.105   #coefficient
        c_z = 1     #coefficient
        c_s = 1     #coefficient
        m = 0.72    #coefficient
        Nu = c*c_z*c_s*(d_ex/pitch)**(-0.54)*(fin_height/pitch)**(-0.14)*Re**m      #Nusselt number
        return Nu
    
    def aNusselt_air(self):
        c = 0.105   #coefficient
        c_z = 1     #coefficient
        c_s = 1     #coefficient
        m = 0.72    #coefficient
        d_ex = self.d_ex
        pitch = self.pitch
        fin_height = self.fin_height
        
        Re = term.Re(term.velFromMas(mass_air,self.hydDiam_air,CP.PropsSI('Dmass','T',T_inside,'P',p_airin,'Air')),self.hydDiam_air,CP.PropsSI('V','T',T_inside,'P',p_airin,'Air'))
        Nu = c*c_z*c_s*(d_ex/pitch)**(-0.54)*(fin_height/pitch)**(-0.14)*Re**m
        return Nu
        
        
        


#convection heat transfer coefficient for air
    def Alfa_air(self,Nusselt,lambda_air,d_ex):
        alfa = (Nusselt*lambda_air)/d_ex
        return alfa
    
#convection heat transfer coefficient for air    
    def aAlfa_air(self):
        alfa = (self.aNusselt_air()*CP.PropsSI('conductivity','T',T_inside,'P',p_airin,'Air'))/self.hydDiam_air
        return alfa
    
    def aFin_module(self):
        a = self.aAlfa_air()
        mod = (2*a/self.lambda_wall/self.delta)
        return mod
    def aFin_efficiency(self):
        mod = self.aFin_module()
        E = m.tanh(mod*self.fin_height)/(mod*self.fin_height)
        return E
    def aFin_surf_eff(self):
        E = self.aFin_efficiency()
        efficiency = E + (1-E)/(self.Area_ext/self.Area_in)
        return efficiency
    
    
#boiling heat transfer of refrigerant inside tube 
    def aAlfa_boil(self,p_boil):
        
        d_i = self.d_in
        A = self.Area_ref
        
        #mass_velocity = 20.0
        
        Q = (1-x_evapin)/2+x_evapin #mean vapor quality along heat exchanger
        sigma = CP.PropsSI('I','P',p_boil,'Q',Q,'R600a')
        rho_gas = CP.PropsSI('DMASS','P',p_boil,'Q',1,'R600a')
        rho_liquid = CP.PropsSI('DMASS','P',p_boil,'Q',0,'R600a')
        mi_liquid = CP.PropsSI('VISCOSITY','P',p_boil,'Q',0,'R600a')
        cp_liquid = CP.PropsSI('CPMASS','P',p_boil,'Q',0,'R600a')
        cp_gas = CP.PropsSI('CPMASS','P',p_boil,'Q',1,'R600a')
        lambda_liquid = CP.PropsSI('CONDUCTIVITY','P',p_boil,'Q',0,'R600a')
        lambda_gas = CP.PropsSI('CONDUCTIVITY','P',p_boil,'Q',1,'R600a')
        mi_gas = CP.PropsSI('VISCOSITY','P',p_boil,'Q',1,'R600a')
        g = 9.81 #earth acceleration
        mass_velocity = mass_flow/A
        
        #loop for heat flux density 
        T_boil = CP.PropsSI('T','P',p_boil,'Q',0,'R600a')
        Twall_in = T_boil+0.0005
        alfa = 100
        convergance = False
        while convergance == False: 
        
            q1 = (T_inside-Twall_in)/(1/self.aAlfa_air()+m.log((self.d_ex/2)/(self.d_in/2))/(2*3.14*self.lambda_wall)) 
        
            alfa_nb = 55*fin_evap.p_r**(0.12-0.4343*m.log(1))*(-m.log(fin_evap.p_r,10))**(-0.55)*fin_evap.M**(-0.5)*q1**(0.67)
            
            epsilon = (Q/rho_gas)*((1+0.12*(1-Q))*(Q/rho_gas+(1-Q)/rho_liquid)+1.18/mass_velocity*(g*sigma*(rho_liquid-rho_gas)/rho_liquid**2)**(0.25)*(1-Q))**(-1)
            
            A_L = A*(1-epsilon)
            theta_strat = 2*(3.14-(m.cos(2*epsilon-1))**(-1))  #not sure how to implement it yet
            
            #theta_strat = 0
            #print (theta_strat)
            #print (p_boil)
            #print (mass_velocity)
            #print (Q)
            #print (epsilon)
            #print (mi_liquid)
            #print (mi_gas)
            #print (x_evapin)
            #print (h_evapin)
            delta = A_L / (d_i/2*(2*3.14-theta_strat))
            
            
            alfa_cb = 0.0133*((4*mass_velocity*(1-Q)*delta)/((1-epsilon)*mi_liquid))**(0.69)*((cp_liquid*mi_liquid)/lambda_liquid)**(0.4)*(lambda_liquid/delta)
            
            alfa_vapor = 0.023*(mass_velocity*Q*d_i/(epsilon*mi_gas))**(0.8)*(cp_gas*mi_gas/lambda_gas)**(0.4)*lambda_gas/d_i
            
            alfa_wet = (alfa_cb**3+alfa_nb**3)**(1/3)
            alfa = (d_i/2*theta_strat*alfa_vapor+d_i/2*(2*3.14-theta_strat)*alfa_wet)/(2*3.14*d_i/2)
            q2 = (Twall_in-T_boil)/(1/alfa)
            eps = abs(q1-q2)
            Twall_in = Twall_in + 0.0005
            #print ('eps = ',eps)
            #print ('q1 = ',q1)
            #print ('q2 = ',q2)
            #print ('alfa = ',alfa)
            #print ('Twall = ',Twall_in)
            if (eps < 20):
                convergance = True 
        alfa = alfa    
        return alfa




    
    def aOHTC(self,p_boil):        #in this kind of solution heat transfer surface is the inner tube surface !!!!
        alfa_p = self.aAlfa_air()
        alfa_r = self.aAlfa_boil(p_boil)
        N = self.aFin_surf_eff()
        OHTC = 1/(self.Area_in/(alfa_p*N*self.Area_ext)+m.log((self.d_ex/2)/(self.d_in/2))/(2*3.14*self.lambda_wall)+1/alfa_r)
        return OHTC
    
    def aOHTC_boil(self,alfa_boil,Area_ext,Area_in):        #in this kind of solution heat transfer surface is the inner tube surface !!!!
        alfa_p = self.aAlfa_air()
        alfa_r = alfa_boil
        N = self.aFin_surf_eff()
        OHTC = 1/(Area_in/(alfa_p*N*Area_ext)+m.log((self.d_ex/2)/(self.d_in/2))/(2*3.14*self.lambda_wall)+1/alfa_r)
        return OHTC
    
    
    
    def NTU(self,p_boil):
        U = self.aOHTC(p_boil)
        A = self.Area_in
        cp = CP.PropsSI('CPMASS','T',T_inside,'P',p_airin,'Air')
        C_min = mass_air*cp
        NTU = U*A/C_min
        return NTU
    
    def NTU_boil(self,A,U):
       
        cp = CP.PropsSI('CPMASS','T',T_inside,'P',p_airin,'Air')
        C_min = mass_air*cp
        NTU = U*A/C_min
        return NTU
    
    def NTU_overh(self,A,U,T_air,T_ref,p_boil):
        cp_ref = CP.PropsSI('CPMASS','T',T_ref,'P',p_boil,'R600a')
        cp_air = CP.PropsSI('CPMASS','T',T_air,'P',p_airin,'Air')
        C1 = cp_ref*mass_flow
        C2 = cp_air*mass_air
        if(C1>C2):
            C_min = C2
        elif(C2>C1):
            C_min = C1
        
        NTU = U*A/C_min
        return NTU
    
    def epsilon_overh(self,NTU,T_ref,T_air,p_boil):
        cp_ref = CP.PropsSI('CPMASS','T',T_ref,'P',p_boil,'R600a')
        cp_air = CP.PropsSI('CPMASS','T',T_air,'P',p_airin,'Air')
        C1 = cp_ref*mass_flow
        C2 = cp_air*mass_air
        if(C1>C2):
            C_min = C2
            C_max = C1
        elif(C2>C1):
            C_min = C1
            C_max = C2
        
        c = C_min/C_max
        epsilon = 1 - m.exp((NTU**0.22)/c*(m.exp(-c*NTU**0.78)-1))
        return epsilon
    
    
    def epsilon_boil(self,NTU,):
        epsilon = 1 - m.exp(-NTU)
        return epsilon
    
    def aHE_capa(self,T_ref,T_air,eps):
        
        Q_max = mass_air*CP.PropsSI('CPMASS','T',T_air,'P',p_airin,'Air')*(T_air- T_ref)
        Q = Q_max * eps
        return Q
    
    def aHE_capaover(self,T_ref,T_air,eps,p_boil):
        cp_ref = CP.PropsSI('CPMASS','T',T_ref,'P',p_boil,'R600a')
        cp_air = CP.PropsSI('CPMASS','T',T_air,'P',p_airin,'Air')
        C1 = cp_ref*mass_flow
        C2 = cp_air*mass_air
        if(C1>C2):
            C_min = C2
        elif(C2>C1):
            C_min = C1
            
        Q = eps*C_min*(T_air-T_ref)
        return Q
    
    
    def aT_airout(self,Q):
        delta_h = Q/mass_air
        h_1 = CP.PropsSI('H','P',p_airin,'T',T_inside,'Air')
        T_airout = CP.PropsSI('T','P',p_airin,'H',h_1-delta_h,'Air')
        return T_airout
    
    def aNusseltoverheat(self,Q_boil,p_boil):
        h_refin = h_evapin + Q_boil/mass_flow
        rho = CP.PropsSI('DMASS','P',p_boil,'H',h_refin,'R600a')
        Prandtl = CP.PropsSI('PRANDTL','P',p_boil,'H',h_refin,'R600a')
        mi = CP.PropsSI('VISCOSITY','P',p_boil,'H',h_refin,'R600a')
        vis = mi/rho
        u = mass_flow/(self.Area_ref*rho)
        Re = term.Re(u,self.d_in,vis)
        Nu = 0.027*Re**(4/5)*Prandtl**(1/3)
        return Nu
    
    def aAlfaoverheat(self,Q_boil,Nu,p_boil):
        h_refin = h_evapin + Q_boil/mass_flow
        alfa = (Nu*CP.PropsSI('conductivity','H',h_refin,'P',p_boil,'R600a'))/self.d_in
        return alfa
    
    
    def aHE_capa_corrected(self):
        n_rows = self.n_rows
        n_columns = self.n_columns
        n_rowso = 0
        p_boili = p_boil
        con = False
        while (con == False):
            #print ('n_rowso= ',n_rowso)
            Area_air_boil = (2*(self.s_1*self.s_2-3.14/4*self.d_ex**2)+3.14*self.d_ex*(self.pitch-self.delta))*(self.tube_length/self.pitch)*n_rows*n_columns
            #print('Aab=',Area_air_boil)
            Area_air_overh = (2*(self.s_1*self.s_2-3.14/4*self.d_ex**2)+3.14*self.d_ex*(self.pitch-self.delta))*(self.tube_length/self.pitch)*n_rowso*n_columns
            #print('Aao=',Area_air_overh)
            Area_ref_boil = 3.14*self.d_in*self.tube_length*n_columns*n_rows
            #print('Arb=',Area_ref_boil)
            Area_ref_overh = 3.14*self.d_in*self.tube_length*n_columns*n_rowso
            #print('Aro=',Area_ref_overh)
            alfa_boil = self.aAlfa_boil(p_boili)
            OHTC_boil = self.aOHTC_boil(alfa_boil,Area_air_boil,Area_ref_boil)
            NTU_boil = self.NTU_boil(Area_ref_boil, OHTC_boil)
            #print ('NRU_boil=',NTU_boil)
            epsilon_boil = self.epsilon_boil(NTU_boil)
            T_boil = CP.PropsSI('T','P',p_boil,'Q',1,'R600a')
            Q_boil = self.aHE_capa(T_boil, T_inside, epsilon_boil)
            #print ('Q_boil=',Q_boil)
            q_boil = Q_boil/mass_flow
            h_refin = h_evapin + q_boil
            T_refin = CP.PropsSI('T','P',p_boil,'H',h_refin,'R600a')
            if n_rowso >0:
                T_airin = self.aT_airout(Q_boil)
                Nusselt_ref = self.aNusseltoverheat(Q_boil,p_boili)
                alfa_ref = self.aAlfaoverheat(Q_boil, Nusselt_ref,p_boili)
                OHTC_over = self.aOHTC_boil(alfa_ref, Area_air_overh, Area_ref_overh)
                NTU_over = self.NTU_overh(Area_ref_overh, OHTC_over, T_airin, T_refin+0.5,p_boili)
                epsilon_over = self.epsilon_overh(NTU_over, T_refin+0.5, T_airin,p_boili)
                Q_over = self.aHE_capaover(T_refin+0.5, T_airin, epsilon_over,p_boili)
            if n_rowso == 0:
                Q_over = 0
            if T_refin > (T_boil+1):
                n_rowso = n_rowso+1
                n_rows = n_rows - 1
            elif T_refin < (T_boil+1):
                con = True
        if (Q_over >0):
            Q = Q_boil + Q_over
        elif (Q_over == 0):
            Q = Q_boil
        return Q

class WoTCondenser: 
    
    '''
    
    This class is based on 'Modeling and optimisation of wire and tube condenser' 
    by P.K. Banasal, T.C. Chin 
    publicated on International Journal of Refrigeration 26 (2003)
    some simplification were used to reduce computational effort
    
    '''
    HEtype = 'Wire on Tube Condenser'

    
    def __init__(self,d_ex,d_in,d_w,pitch_tube,pitch_wire,tube_length,n_columns):
        self.d_ex = d_ex*10**(-3)   
        self.d_in = d_in*10**(-3)                
        self.pitch_tube = pitch_tube*10**(-3)
        self.pitch_wire = pitch_wire*10**(-3)          
        self.tube_length = tube_length*10**(-3)
        self.n_columns = n_columns
        
        self.Area_ref = 3.14*self.d_in**2/4
        self.lambda_wall = 45
        self.d_w = d_w*10**(-3)
        self.Area_ext = (3.14*self.d_ex*(self.tube_length +3.14*self.pitch_tube)*self.n_columns +3.14*self.d_w*self.pitch_tube*(self.tube_length/self.pitch_wire)*2)*2
        self.Area_in = 3.14*self.d_in*(self.tube_length+3.14*self.pitch_tube/2)*self.n_columns*2
        self.s_t = (self.pitch_tube - self.d_ex)/self.d_ex
        self.s_w = (self.pitch_wire-self.d_w)/self.d_w
        self.H = self.n_columns*self.pitch_tube
        
        
    def aRayleigh(self,T_tube):
        T_inf = T_env
        beta = 1/T_inf
        ro = CP.PropsSI('DMASS','P',101300,'T',T_inf,'Air')
        cp = CP.PropsSI('CP0MASS','P',101300,'T',T_inf,'Air')
        mi = CP.PropsSI('VISCOSITY','P',101300,'T',T_inf,'Air')
        lambda_air = CP.PropsSI('CONDUCTIVITY','P',101300,'T',T_inf,'Air')
        g = 9.81
        Ra = (beta*ro**2*cp/(mi*lambda_air))*g*(T_tube-T_inf)*self.H**3
        #print ('Ra= ',Ra)
        return Ra
    
    def aPhi(self,T_tube):
        T_inf = T_env
        Phi = (28.2/self.H)**(0.4)*self.s_w**(0.9)*self.s_t**(-1)+(28.2/self.H)**(0.8)*(264/(T_tube-T_inf))**(0.5)*self.s_w**(-1.5)*self.s_t**(-0.5)
        #print ('Phi=',  Phi)
        return Phi
   
    def aNusselt(self,T_tube):
        Nu = 0.66*(self.aRayleigh(T_tube)*self.H/self.d_ex)**(0.25)*(1-(1-0.45*(self.d_ex/self.H)**(0.25))*m.exp(-self.s_w/self.aPhi(T_tube)))
        #print ('Nu=', Nu)
        return Nu
    
    def aAlfaAir(self,T_tube):
        T_inf = T_env
        lambda_air = CP.PropsSI('CONDUCTIVITY','P',101300,'T',T_inf,'Air')
        alfa = self.aNusselt(T_tube)*lambda_air/self.H
        alfa = alfa
        return alfa
    
    def aFin_module(self,alfa_w):
        module = (4*alfa_w/(self.lambda_wall*self.d_w))**(1/2)
        return module
    
    def aFin_efficiency(self,alfa_w):
        module = self.aFin_module(alfa_w)
        E = m.tanh(module*self.pitch_tube/2)/(module*self.pitch_tube/2)
        return E
   
    def aT_ex(self,T_tube,T_w,alfa_w):
        T_inf = T_env
        GP = 2*(self.pitch_tube/self.d_ex)*(self.d_w/self.pitch_wire)
        T_ex = (T_tube+GP*self.aFin_efficiency(alfa_w)*(T_tube-T_inf)+GP*T_inf)/(1+GP)
        return T_ex
    
    def aAlfaAir_vertical(self,T_tube):
        T_inf = T_env
        alfa = 0.27*((T_inf-T_tube)/self.d_ex)**(0.25)
        return alfa
    
    def aMassVel(self): #this will be used to identify multiphase flow pattern sometime in future
        G = mass_flow/self.Area_ref
        x = 0.5
        p = p_cond
        rho_gas = CP.PropsSI('DMASS','P',p,'Q',1,'R600a')
        rho_liq = CP.PropsSI('DMASS','P',p,'Q',0,'R600a')
        j_g = x*G/(9.81*self.d_in*rho_liq*(rho_liq-rho_gas))
        return j_g
    
    def aLMparameter(self):     #Lockhart Martinelli parameter will be used to identify flow patern
        p = p_cond
        x = 0.5
        rho_gas = CP.PropsSI('DMASS','P',p,'Q',1,'R600a')
        rho_liq = CP.PropsSI('DMASS','P',p,'Q',0,'R600a')
        mi_liquid = CP.PropsSI('VISCOSITY','P',p,'Q',0,'R600a')
        mi_gas = CP.PropsSI('VISCOSITY','P',p,'Q',1,'R600a')
        X_u = ((1-x)/x)**(0.875)*(rho_liq/rho_gas)**(0.5)*(mi_liquid/mi_gas)**(0.5)
        return X_u
    
    def aAlfaCond(self): #Condensation for annular flow in horizontal tubes, for now it is fixed pattern
        x = 0.5
        p = p_cond
        G = mass_flow/self.Area_ref
        lambda_liquid = CP.PropsSI('CONDUCTIVITY','P',p,'Q',0,'R600a')
        Prandtl = CP.PropsSI('PRANDTL','P',p,'Q',0,'R600a')
        rho_gas = CP.PropsSI('DMASS','P',p,'Q',1,'R600a')
        rho_liq = CP.PropsSI('DMASS','P',p,'Q',0,'R600a')
        mi_liquid = CP.PropsSI('VISCOSITY','P',p,'Q',0,'R600a')
        mi_gas = CP.PropsSI('VISCOSITY','P',p,'Q',1,'R600a')
        Re_v = G*x*self.d_in/mi_gas
        Re_l = G*(1-x)*self.d_in/mi_liquid
        Re_eq = Re_v*(mi_gas/mi_liquid)*(rho_liq/rho_gas)**(0.5)+Re_l
        Nu = 0.05*Re_eq**(0.8)*Prandtl**(0.3)
        alfa = Nu * lambda_liquid/self.d_in
        
        return alfa
    
    def aOHTC(self,alfa_air,alfa_cond,eta):
        OHTC = 1/(self.Area_in/(alfa_air*eta*self.Area_ext)+m.log((self.d_ex/2)/(self.d_in/2))/(2*3.14*self.lambda_wall)+1/alfa_cond)
        return OHTC
    
    
    def aAlfaRad(self, T_ex):
        sigma = 5.67*10**(-8)
        epsilon = 0.91
        alfa = sigma * epsilon *(T_ex**4-T_env**4)/(T_ex - T_env)
        return alfa
    
    
    def aHE_capa(self):  #main loop for evaluation condenser capacity 
        
        T_cond = CP.PropsSI('T','P',p_cond,'Q',0,'R600a')
        T_tube = T_cond - 0.01
        alfa_w = 1
        conv = False
        convergance = False
        while (convergance == False):
            
            while (conv == False):
                    #print('alfa_w=',alfa_w)
                    T_inf = T_env
                    eta = self.aFin_efficiency(alfa_w)
                    T_w = eta*(T_tube - T_inf) + T_inf
                    T_ex = self.aT_ex(T_tube, T_w, alfa_w)
                    alfa_air = self.aAlfaAir(T_ex)
                    alfa_rad = self.aAlfaRad(T_ex)
                    alfa_all = alfa_air + alfa_rad
                    convi = abs(alfa_all - alfa_w)
                    #print(convi)
                    if abs(convi) > 0.1:
                        alfa_w = alfa_w + 0.01
                    if (abs(convi) < 0.1) & (convi >0.01):
                        alfa_w = alfa_all
                    
                    elif abs(convi) < 0.01:
                        conv = True   
                       
            alfa_cond = self.aAlfaCond()
            #OHTC = self.aOHTC(alfa_air, alfa_cond, eta)
            A = self.Area_ext
            #Q = OHTC*A*(T_cond-T_inf)
            R = 1/(alfa_cond*self.Area_in)+m.log(self.d_ex/2/self.d_in/2)/(2*3.14*self.lambda_wall*self.tube_length*self.n_columns)+1/(self.Area_ext*alfa_air)
            #R = 1/(1/alfa_air/self.Area_ext-0.005/self.lambda_wall+1/alfa_cond/self.Area_in)
            #print ('R = ',R)
            #Q = OHTC*A*(T_cond-T_inf)
            Q = (T_cond-T_inf)/R
            #print ('OHTCcond=',OHTC)
            #print ('Q=',Q)
            T = T_cond - Q*(1/(alfa_cond*self.Area_in)+m.log(self.d_ex/2/self.d_in/2)/(2*3.14*self.lambda_wall*self.tube_length*self.n_columns))
            #T = T_cond - Q/R/A
            #print ('T=',T)
            #print ('Tcond=',T_cond)
            #print ('Ttube=',T_tube)
            #print ('alfa_air= ',alfa_air)
            #print ('alfa_cond=', alfa_cond)
            delta_T = T_tube - T
            #print ('deltaT= ',delta_T)
            if abs(delta_T) >= 0.01:
                T_tube = T_tube - 0.01
                conv = False
            elif abs(delta_T) <0.01:
                convergance = True
                
        return Q

class compressor:
    
    speed = 2400    #actual compressor speed
    
    def __init__(self,N,V,start_loss,max_eff):
        self.N = N                          #rotational speed, rpm
        self.V = V                          #volumetric displacement, cm3
        self.start_loss = start_loss        #start loss, W
        self.max_eff = max_eff              #max efficency for working point , -
        #self.p_boil = p_boil
        
    
    
        
    def aEfficiency(self):
        x = compressor.speed
        y = self.max_eff
        
        if (compressor.speed == self.N):
            eff = y
        elif (compressor.speed > self.N):
            b = 1
            a = (y-b)/self.N
            eff = a*x + b
        elif (compressor.speed < self.N):
            b = 0.3
            a = (y-b)/self.N
            eff = a*x + b
            
        return eff
  
    def aMassFlow(self,p_boil):
        ro = CP.PropsSI('DMASS','P',p_boil,'T',T_inlet,'R600a')
        volFlow = self.V*compressor.speed/60*0.000001
        massFlow = volFlow*ro
        return massFlow
    
    def aWork(self,p_boil,p_cond,T_inlet):
        s = CP.PropsSI('S','T',T_inlet,'P',p_boil,'R600a')
        h_inlet = CP.PropsSI('H','T',T_inlet,'P',p_boil,'R600a')
        h_outlet_perf = CP.PropsSI('H','S',s,'P',p_cond,'R600a')
        h_outlet = h_inlet + (h_outlet_perf-h_inlet)/self.aEfficiency()
        work = h_outlet - h_inlet
        return work
    
    def aPower(self,p_boil,p_cond,T_inlet):
        Power = self.aWork(p_boil,p_cond,T_inlet)*self.aMassFlow(p_boil)
        return Power
    
    def aT_outlet(self,p_boil,p_cond,T_inlet):
        h_inlet = CP.PropsSI('H','T',T_inlet,'P',p_boil,'R600a')
        deltaH = self.aWork(p_boil,p_cond,T_inlet)
        T = CP.PropsSI('T','H',h_inlet+deltaH,'P',p_cond,'R600a')
        return T

class Wall:
    
    def __init__(self,d1,d2,d3,lambda1,lambda2,lambda3,h,w):
        self.d1 = d1*10**(-3)    #wall thickness 
        self.d2 = d2*10**(-3) 
        self.d3 = d3*10**(-3) 
        self.lambda1 = lambda1  #wall conductivity
        self.lambda2 = lambda2
        self.lambda3 = lambda3
        self.h = h    #height of wall
        self.w = w      #width of wall
        
    def aAlfaOut(self):
        d = self.h
        T = T_env
        beta = 1/T
        Prandtl = CP.PropsSI('Prandtl','T',T,'P',101300,'Air')
        density = CP.PropsSI('DMASS','T',T,'P',101300,'Air')
        visco_kin = CP.PropsSI('V','T',T,'P',101300,'Air')/density 
        grashof = 9.81*d**3*beta*5/visco_kin**2
        lambda_air = CP.PropsSI('CONDUCTIVITY','T',T,'P',101300,'Air')
        
        if (grashof*Prandtl<= 10**(-3)):
            C = 0.5
            n = 0
        elif (grashof*Prandtl> 10**(-3) and grashof*Prandtl<= 500):
            C = 1.18
            n = 0.125
        elif (grashof*Prandtl> 500 and grashof*Prandtl<= 2*10**7):
            C = 0.54
            n = 0.25
        elif (grashof*Prandtl> 2*10**7 and grashof*Prandtl<= 1*10**13):
            C = 0.135
            n = 1/3
        
        Nu = C*(grashof*Prandtl)**n
        alfa = Nu * lambda_air / d
        return alfa
    
    def aAlfaIn(self):
        d = self.h
        T = T_inside
        beta = 1/T
        Prandtl = CP.PropsSI('Prandtl','T',T,'P',101300,'Air')
        density = CP.PropsSI('DMASS','T',T,'P',101300,'Air')
        visco_kin = CP.PropsSI('V','T',T,'P',101300,'Air')/density 
        grashof = 9.81*d**3*beta*5/visco_kin**2
        lambda_air = CP.PropsSI('CONDUCTIVITY','T',T,'P',101300,'Air')
        
        if (grashof*Prandtl<= 10**(-3)):
            C = 0.5
            n = 0
        elif (grashof*Prandtl> 10**(-3) and grashof*Prandtl<= 500):
            C = 1.18
            n = 0.125
        elif (grashof*Prandtl> 500 and grashof*Prandtl<= 2*10**7):
            C = 0.54
            n = 0.25
        elif (grashof*Prandtl> 2*10**7 and grashof*Prandtl<= 1*10**13):
            C = 0.135
            n = 1/3
        
        Nu = C*(grashof*Prandtl)**n
        alfa = Nu * lambda_air / d
        return alfa
    
    def aOHTC(self):
        
        alfa1 = self.aAlfaIn()
        alfa2 = self.aAlfaOut()
        OHTC = 1/(1/alfa1 + self.d1/self.lambda1 + + self.d2/self.lambda2 + + self.d3/self.lambda3 + 1/alfa2)
        return OHTC
    def aHeat(self):
        k = self.aOHTC()
        
        Q = self.w*self.h*k*(T_env-T_inside)
        return Q

class capillary:  #not working atm
    
    def __init__(self,D,L,instances,Q):
        self.D = D*0.001    #converts diameter to SI
        self.L = L*0.001    #converts length to SI
        #self.d = self.D
        #self.l = self.L/capillary.num_instances
        self.instances = instances
        self.l = self.L/self.instances
        self.Q = Q
        
    def aPressDropSimple(self):
        ro = ro_inlet_capi
        mflow = mass_flow
        p = p_cond
        i = 1
        dQ = self.Q/self.instances
        q = dQ/mass_flow
        h = h_condout
        while i < self.instances:
            #print (p)
            h = h - q
            v = term.velFromMas(mflow,self.D,ro)
            mi = CP.PropsSI('VISCOSITY','P',p,'H',h,'R600a')/ro
            Re = term.Re(v,self.D,mi) 
            lambdapd = term.pdCoeff(Re)
            pd = term.pressDrop(lambdapd,v,self.D,ro)*self.l
            #print (pd)
            p = p-pd
            ro = CP.PropsSI('DMASS','P',p,'H',h,'R600a')
            
            #print (p) 
            i = i + 1
            #print (i,'/',self.instances)
            #print (Re)
            #print ('h=', h)
            #print ('h_inlet=',h_condout)
            
        #Q = (h_condout-h)*mass_flow   
        #x = CP.PropsSI('Q','P',p,'H',h,'R600a')
        return p

class GlassPlate:
    
    def __init__(self,h,w,l):
        self.h = h*10**(-3) #height of plate
        self.w = w*10**(-3) #width of plate
        self.l = l*10**(-3) #length of plate
        
    def aMassCP(self):
        rho = 2500
        cp = 840
        V = self.h*self.w*self.l
        m = V*rho
        MassCP = m*cp
        return MassCP

class PlasticBox:
    
    def __init__(self,h,w,l,t):
        self.h = h*10**(-3) #height of plate
        self.w = w*10**(-3) #width of plate
        self.l = l*10**(-3) #length of plate
        self.t = t*10**(-3) #wall thickness
        
    def aMassCP(self):
        rho = 920
        cp = 1.67
        V = self.t*self.w*self.l + 2*self.h*self.w*self.t + 2*self.h*self.l*self.t
        m = V*rho
        MassCP = m*cp
        return MassCP




'''
Construction of all freezer components
'''
height = 1.6                #height of compartment (inside)
width = 0.55                 #width of compartment (inside)
deep = 0.55                  #deep of compartment (inside)
Q_capi = 10
refrigerant = 'R600a'
compressor = compressor(2400,7,100,0.8)
evaporator = fin_evap(6,5.2,3.5,450,16,2,22,22)
condenser = WoTCondenser(5, 4, 1, 50, 4, 600, 27)
capillary = capillary(0.6,2500,1000,Q_capi)
door = Wall(0.6,120,4,40,0.03,2,1.6,0.5)
leftWall = Wall(0.4,100,4,40,0.03,2,1.6,0.55)
rightWall = Wall(0.4,100,4,40,0.03,2,1.6,0.55)
backWall = Wall(0.4,100,4,40,0.03,2,1.6,0.55)
topWall = Wall(0.4,80,4,40,0.03,2,0.55,0.55)
bottomWall = Wall(0.4,80,4,40,0.03,2,0.55,0.55)
box = PlasticBox(200, 550, 550, 5)
glass = GlassPlate(10, 550, 550)
V_compartment = height*width*deep           #Volume of compartment
p_boil = 0.59*10**5

'''
Creation of output files
'''
presscondoutput = open("HighPressure.txt","w+")
presscondoutput.write("This file contains condensation pressure for every time-step \n")
presscondoutput.close()

pressevapoutput = open("LowPressure.txt","w+")
pressevapoutput.write("This file contains evaporation pressure for every time-step \n")
pressevapoutput.close()

powerconsumoutput = open("PowerConsumption.txt","w+")
powerconsumoutput.write("This file contains power consumption for every time-step \n")
powerconsumoutput.close()

coolingcapaoutput = open("CoolingCapacity.txt","w+")
coolingcapaoutput.write("This file contains cooling capacity for every time-step \n")
coolingcapaoutput.close()

heatinleakoutput = open("Heatinleak.txt","w+")
heatinleakoutput.write("This file contains heat inleaks for every time-step \n")
heatinleakoutput.close()

tempoutput = open("Temperature.txt","w+")
tempoutput.write("This file contains freezer's temperature for every time-step \n")
tempoutput.close()

timeoutput = open("Time.txt","w+")
timeoutput.write("This file contains time values \n")
timeoutput.close()

'''
Start of the main time loop
'''
rest = False

while (time < end_time):
    time_off = 0
    
    #print(rest)
    #print ('T_insC= ',T_inside-273)
    if (T_inside > T_var1 and rest == False):
        
        '''
        Begin of refrigeration cycle capacity
        '''
        
        x_conv = False
        while (x_conv == False): #Condenser Pressure Loop
        
            mass_flow = compressor.aMassFlow(p_boil)
            #print('massflow=',mass_flow)
            T_compout = compressor.aT_outlet(p_boil,p_cond,T_inlet)
            power_usage = compressor.aPower(p_boil,p_cond,T_inlet)
            '''
            Here should be used frame heater but for now it does not exist
            '''
            h_condin = CP.PropsSI('H','P',p_cond,'Q',1,refrigerant)
            Q_cond = condenser.aHE_capa()
            #print('Qcond = ',Q_cond)
            q_cond = Q_cond/mass_flow                               #condenser capacity J/kg 
            h_condout = h_condin - q_cond
            x_condout = CP.PropsSI('Q','P',p_cond,'H',h_condout,refrigerant)
            #print('xcond=',x_condout)
            #print('p_cond=',p_cond)
            if (x_condout == (-1) ):
                x_conv = True
            elif (x_condout != (-1)):
                p_cond = p_cond + 0.5*10**4
        
        '''
        Capillary tube will be added someday for now it is not working
        '''
        p_boil = 0.59*10**5
        p_conv = False
        while (p_conv == False): #Evaporator Pressure Loop
            #print ('P_boil=',p_boil)
            #ro_inlet_capi = CP.PropsSI('DMASS','P',p_cond,'H',h_condout,refrigerant)
            #p_boil = capillary.aPressDropSimple()
            #print('pboil = ',p_boil)
            q_capi = Q_capi/mass_flow
            h_evapin = h_condout-q_capi
            x_evapin = CP.PropsSI('Q','P',p_boil,'H',h_evapin,refrigerant)
            #print('x_evapin=',x_evapin)
            if x_evapin == -1:
                x_evapin = 0
                h_evapin = CP.PropsSI('H','P',p_boil,'Q',x_evapin,refrigerant)
            
            Q_evap = evaporator.aHE_capa_corrected()
            #print('Q_evap',Q_evap)
            q_evap = Q_evap/mass_flow
            h_evapout = h_evapin + q_evap
            x_evapout = CP.PropsSI('Q','P',p_boil,'H',h_evapout,refrigerant)
            #print('x_evapout=',x_evapout)
            q_suctionline = q_capi
            h_compinlet = h_evapout + q_suctionline
            T_inlet = CP.PropsSI('T','P',p_boil,'H',h_compinlet,refrigerant)
            #print('Tinlet',T_inlet)
            h_x1 = CP.PropsSI('H','P',p_boil,'Q',1,refrigerant)
            if (h_compinlet > (h_x1+200)):
                p_conv = True
            elif (h_compinlet < (h_x1+200)):
                p_boil = p_boil - 0.005*10**5
            
            
        '''
        End of refrigeration cycle calculations 
        '''
        '''
        Heat Balance of compartment
        '''
        
        Q_rw = rightWall.aHeat()
        Q_lw = leftWall.aHeat()
        Q_door = door.aHeat()
        Q_bw = backWall.aHeat()
        Q_tw = topWall.aHeat()
        Q_botw = bottomWall.aHeat()
        Q_inleak = Q_rw + Q_lw + Q_door + Q_bw + Q_tw + Q_botw
        #print('Qin=',Q_inleak)
        Q_real = Q_inleak - Q_evap
        #print('dQ=',Q_real)
        rho_air = CP.PropsSI('DMASS','P',p_airin,'T',T_inside,'Air')
        mass_air = V_compartment*rho_air
        cp_air = CP.PropsSI('CPMASS','P',p_airin,'T',T_inside,'Air')
        dT = Q_real/(mass_air*cp_air + 5*box.aMassCP() + 7*glass.aMassCP())
        #print ('dT = ',dT)
        T_inside = T_inside + dT
        #print('T_inside=',T_inside)
        print('time = ',time,'/',end_time)
        time = time + delta_time
        
        if (T_inside < T_var1 and rest == False):
            rest = True
            
    elif (T_inside < T_var2 and rest == True):
        
        Q_evap = 0
        Q_rw = rightWall.aHeat()
        Q_lw = leftWall.aHeat()
        Q_door = door.aHeat()
        Q_bw = backWall.aHeat()
        Q_tw = topWall.aHeat()
        Q_botw = bottomWall.aHeat()
        Q_inleak = Q_rw + Q_lw + Q_door + Q_bw + Q_tw + Q_botw
        #print('Qin=',Q_inleak)
        Q_real = Q_inleak
        #print('dQ=',Q_real)
        rho_air = CP.PropsSI('DMASS','P',p_airin,'T',T_inside,'Air')
        mass_air = V_compartment*rho_air
        cp_air = CP.PropsSI('CPMASS','P',p_airin,'T',T_inside,'Air')
        dT = Q_real/(mass_air*cp_air + 5*box.aMassCP() + 7*glass.aMassCP())
        T_inside = T_inside + dT
        #print('T_inside=',T_inside)
        print('time = ',time,'/',end_time)
        time = time + delta_time
        
        if (T_inside > T_var2 and rest == True ):
            rest = False
        
    '''
    Writing output values
    '''
    if (time%30 == 0 or time == 1):
        
        p_cond_str = str(p_cond)
        presscondoutput = open("HighPressure.txt","a")
        presscondoutput.write(p_cond_str)
        presscondoutput.write("\n")
        presscondoutput.close()
        
        p_boil_str = str(p_boil)
        pressevapoutput = open("LowPressure.txt","a")
        pressevapoutput.write(p_boil_str)
        pressevapoutput.write("\n")
        pressevapoutput.close()
        
        power_usage_str = str(power_usage)
        powerconsumoutput = open("PowerConsumption.txt","a")
        powerconsumoutput.write(power_usage_str)
        powerconsumoutput.write("\n")
        powerconsumoutput.close()
        
        Q_evap_str = str(Q_evap)
        coolingcapaoutput = open("CoolingCapacity.txt","a")
        coolingcapaoutput.write(Q_evap_str)
        coolingcapaoutput.write("\n")
        coolingcapaoutput.close()
        
        Q_inleak_str = str(Q_inleak)
        heatinleakoutput = open("Heatinleak.txt","a")
        heatinleakoutput.write(Q_inleak_str)
        heatinleakoutput.write("\n")
        heatinleakoutput.close()
        
        T_inside_str = str(T_inside)
        tempoutput = open("Temperature.txt","a")
        tempoutput.write(T_inside_str)
        tempoutput.write("\n")
        tempoutput.close()

        time_str = str(time)
        timeoutput = open("Time.txt","a")
        timeoutput.write(time_str)
        timeoutput.write("\n")
        timeoutput.close()







    
    
    
    