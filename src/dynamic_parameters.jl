# module Basic  

export Bcr4bp_Aux  

"""
    Bcr4bp_Aux() -> NamedTuple

返回地月—太阳双圆四体模型的量纲参数、无量纲参数和归一化尺度。该返回值也是本包
动力学函数与事件回调的参数 `p`。
"""
function Bcr4bp_Aux()
        # ------------------------- Basic Constants -------------------------  

    
        # Earth gravitational parameter  
        mu_E = 398600.435507 
    
        # Moon gravitational parameter  
        mu_M =  4902.800118 
    
        # Sun gravitational parameter  
        mu_S = 132712440041.279419 
    
        # Combined masses  
        mu_EM = mu_E + mu_M     # Earth-Moon system mass [kg]  
        mu_SB1 = mu_E + mu_M + mu_S  # S and SB1 system mass [kg]  

        ## ------------------------- EMRot -------------------------  
        # Earth-Moon distance [km]  
        l_EM = 384405.  
    
        # Earth-Moon time unit [s]  
        t_EM = sqrt(l_EM^3 / (mu_EM))  
    
        # Earth-Moon velocity unit [km/s]  
        v_EM = l_EM / t_EM  
    
        ## ------------------------- SB1Rot -------------------------  
        # Sun-SB1 distance [km]  
        l_SB1 = 149597870.7
    
        # Sun-SB1 time unit [s]  
        t_SB1 = sqrt(l_SB1^3 / (mu_SB1))  
    
        # Sun-SB1 velocity unit [km/s]  
        v_SB1 = l_SB1 / t_SB1  


        r_p1    = 6378.0 
        r_p2    = 1937.4

        dim = (  
            mu_E    = mu_E,  
            mu_M    = mu_M,  
            mu_S    = mu_S,  
            mu_EM    = mu_EM,  
            mu_SB1   = mu_SB1,  
            r_p1    = r_p1 ,
            r_p2    = r_p2 ,
            EMRot_l = l_EM,
            EMRot_v = v_EM,
            EMRot_t = t_EM,
            SB1Rot_l = l_SB1,
            SB1Rot_v = v_SB1,  
            SB1Rot_t = t_SB1,  
        )  
    
        ## ------------------------- EMRot Normalization -------------------------  

        μ =  mu_M / (mu_EM)
        d2lim = min(0.9, 6(μ/3)^(1/3))

        # d2lim = 5

        ##

        EMRot = (  
            μ         = μ,  
            mus        = mu_S / mu_EM,  
            aem        = 1,  
            as         = l_SB1 / l_EM,  
            ws         = sqrt((1 + mu_S / mu_EM) / (l_SB1 / l_EM)^3) - 1,  
            rMag_b1_e  = mu_M / (mu_EM),  
            rMag_b1_m  = 1 - (mu_M / (mu_EM)),  
            rMag_b1_s  = l_SB1 / l_EM,  
            rMag_b2_b1 = (l_SB1 / l_EM) * (mu_S / mu_EM) / (mu_S / mu_EM + 1),  
            rMag_b2_s  = (l_SB1 / l_EM) * 1 / (mu_S / mu_EM + 1)  , 
            d2lim = d2lim,
            r_p2 = r_p2 / l_EM,
            r_p1 = r_p1 / l_EM,

        )  
    
        ## ------------------------- SB1Rot Normalization -------------------------  
        SB1Rot = (  
            μ         = (mu_EM) / (mu_SB1),  
            mus        = mu_S / (mu_SB1),  
            aem        = l_EM / l_SB1,  
            as         = 1,  
            wm         = 1 / sqrt((mu_S / mu_EM + 1) / (l_SB1 / l_EM)^3) - 1,  
            rMag_b2_s  = 1 / (mu_S / mu_EM + 1),  
            rMag_b2_b1 = (mu_S / mu_EM) / (mu_S / mu_EM + 1),  
            rMag_b1_e  = 1 / (l_SB1 / l_EM) * (mu_M / (mu_EM)),  
            rMag_b1_m  = 1 / (l_SB1 / l_EM) * (1 - mu_M / (mu_EM))  
        )  

    
        ## ------------------------- Save Data into aux -------------------------  
        aux = (dim = dim ,
        EMRot = EMRot,
        SB1Rot = SB1Rot
        ) 

        return aux  
    end


    

# end # module
