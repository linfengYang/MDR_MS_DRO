Fully Adaptive Distributionally Robust Multi-stage Framework for Uncertain Unit Commitment Based on Mixed Decision Rules

================

What about this project/study?

      
      With growing penetration of wind power into the power grid, while achieving low cost sustainable electricity supply, it also 
introduces technical challenges with the associated intermittency. This paper proposes a fully adaptive Wasserstein-based 
distributionally robust multi-stage frame-work based on mixed decision rules (MDR) for uncertain unit commitment problem (UUC) to 
better adapt wind pow-er respecting non-anticipativity both in unit state decision and dispatch process. Comparing with the existing 
multi-stage model, the proposed framework introduces an im-proved MDR to handle all decision variables to expand feasible region, thus,
this framework can obtain various typical models by adjusting the number of relevant periods of decision variables. As a result, our 
model can find feasi-ble solution to some problems that are not feasible in the traditional models while finding better solution to
feasible problems. The proposed model is reformulated with ad-vanced optimization method and improved MDR to form the mixed integer 
linear programming (MILP) model to address the computational intractability. The effectiveness and efficiency of the proposed model
have been validated with case studies using IEEE benchmark systems.


User Guide
-----------

The description of implement code files 

    BDR_UCC_DRO : Our proposed DR_MDR framework.  
    
    TS_UCC_DRO : The typical multi-stage DRO model.
   
    ReadDataSCUC :  Read the data.
    
    ReadWindData :  Read the wind data.
    
    SCUC_nodeY :  Construct network admittance matrix.
    
    L_hat :  Continuous lifting operator.
    
    L_tine :  Binary lifting operator.
    
    Sumt :  Calculate number of periods.
    
    U_matrix :  Construct vertices matrix.
    
    Zbin :  Construct binary coefficient vector.
    
    Zboth :  Construct mixed coefficient vector.
    
    Zhat :  Construct continuous coefficient vector.
    
    Zhat_nonb :  Construct continuous coefficient vector for typical models.
    
    vers :  Construct vertices set.
    
    vers_nonb :  Construct vertices set for typical models.
    
    data: All IEEE datas in the simulations.
    
    SCUC_X_Y ：IEEE X-bus Y-periods test system.
    
    Wind_power : Historal wind data.



Prerequisite:
-----------

    Matlab R2018a
    Cplex 12.7.1
    Mosek 9.2




Publication:
-----------
    If you use our study in academic work then please consider citing our papers.




About Us 
-----------
    Authors：Lingfeng Yang (ylf@gxu.edu.cn),Ying Yang (907803678@qq.com),Zhaoyang Dong
    Team：www.scholat.com/team/eidp
    Webpage: http://jians.gxu.edu.cn/default.do
