Help on package quanestimation:

NAME
    quanestimation - Top-level package for quanestimation.

PACKAGE CONTENTS
    Adaptive (package)
    AsymptoticBound (package)
    BayesianBound (package)
    Common (package)
    ComprehensiveOpt (package)
    ControlOpt (package)
    Dynamics (package)
    MeasurementOpt (package)
    Resources (package)
    StateOpt (package)

CLASSES
    builtins.object
        quanestimation.Adaptive.adaptMZI.adaptMZI
        quanestimation.Adaptive.adaptive.adaptive
        quanestimation.Dynamics.dynamics.Lindblad
    quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem(builtins.object)
        quanestimation.ComprehensiveOpt.AD_Compopt.AD_Compopt
        quanestimation.ComprehensiveOpt.DE_Compopt.DE_Compopt
        quanestimation.ComprehensiveOpt.PSO_Compopt.PSO_Compopt
    quanestimation.ControlOpt.ControlStruct.ControlSystem(builtins.object)
        quanestimation.ControlOpt.DDPG_Copt.DDPG_Copt
        quanestimation.ControlOpt.DE_Copt.DE_Copt
        quanestimation.ControlOpt.GRAPE_Copt.GRAPE_Copt
        quanestimation.ControlOpt.PSO_Copt.PSO_Copt
    quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem(builtins.object)
        quanestimation.MeasurementOpt.AD_Mopt.AD_Mopt
        quanestimation.MeasurementOpt.DE_Mopt.DE_Mopt
        quanestimation.MeasurementOpt.PSO_Mopt.PSO_Mopt
    quanestimation.StateOpt.StateStruct.StateSystem(builtins.object)
        quanestimation.StateOpt.AD_Sopt.AD_Sopt
        quanestimation.StateOpt.DDPG_Sopt.DDPG_Sopt
        quanestimation.StateOpt.DE_Sopt.DE_Sopt
        quanestimation.StateOpt.NM_Sopt.NM_Sopt
        quanestimation.StateOpt.PSO_Sopt.PSO_Sopt
    
    class AD_Compopt(quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem)
     |  AD_Compopt(savefile=False, Adam=False, psi0=[], ctrl0=[], measurement0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, seed=1234, eps=1e-08)
     |  
     |  Method resolution order:
     |      AD_Compopt
     |      quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  SC(self, W=[], M=[], target='QFIM', LDtype='SLD')
     |      Description: use DE algorithm to optimize states and control coefficients.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, Adam=False, psi0=[], ctrl0=[], measurement0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, seed=1234, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      
     |      psi0:
     |         --description: initial guesses of states (kets).
     |         --type: array
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      measurement0:
     |         --description: a set of POVMs.
     |         --type: list (of vector)
     |      
     |      savefile:
     |          --description: True: save the states (or controls, measurements) and the value of the
     |                               target function for each episode.
     |                         False: save the states (or controls, measurements) and all the value
     |                                 of the target function for the last episode.
     |          --type: bool
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem:
     |  
     |  CM(self, rho0, W=[])
     |      Description: use DE algorithm to optimize control coefficients and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      rho0:
     |          --description: initial state.
     |          --type: density matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SCM(self, W=[])
     |      Description: use DE algorithm to optimize states, control coefficients and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SM(self, W=[])
     |      Description: use DE algorithm to optimize states and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[], ctrl_bound=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      psi0:
     |         --description: initial state.
     |         --type: vector
     |      
     |      measurement0:
     |         --description: a set of POVMs.
     |         --type: list (of vector)
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save_meas(self)
     |  
     |  load_save_state(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class AD_Mopt(quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem)
     |  AD_Mopt(mtype, minput, savefile=False, Adam=False, measurement0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      AD_Mopt
     |      quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, W=[])
     |      Description: use differential evolution algorithm to update the measurements that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, mtype, minput, savefile=False, Adam=False, measurement0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, seed=1234, load=False, eps=1e-08)
     |       ----------
     |       Inputs
     |       ----------
     |      savefile:
     |           --description: True: save the measurements and the value of the target function
     |                                for each episode.
     |                          False: save the measurements and all the value of the target
     |                                 function for the last episode.
     |           --type: bool
     |      
     |       measurement0:
     |          --description: a set of POVMs.
     |          --type: list (of vector)
     |      
     |       eps:
     |           --description: calculation eps.
     |           --type: float
     |      
     |       notes: the Weyl-Heisenberg covariant SIC-POVM fiducial state of dimension $d$
     |              are download from http://www.physics.umb.edu/Research/QBism/solutions.html.
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, rho0, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class AD_Sopt(quanestimation.StateOpt.StateStruct.StateSystem)
     |  AD_Sopt(savefile=False, Adam=False, psi0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      AD_Sopt
     |      quanestimation.StateOpt.StateStruct.StateSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, Adam=False, psi0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, seed=1234, load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      savefile:
     |          --description: True: save the states and the value of target function for each episode .
     |                         False: save the states and the value of target function for the last episode.
     |          --type: bool
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      seed:
     |          --description: random seed.
     |          --type: int
     |      
     |      eps:
     |          --description: machine eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class DDPG_Copt(quanestimation.ControlOpt.ControlStruct.ControlSystem)
     |  DDPG_Copt(savefile=False, max_episode=500, layer_num=3, layer_dim=200, seed=1234, ctrl0=[], load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      DDPG_Copt
     |      quanestimation.ControlOpt.ControlStruct.ControlSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   HCRB.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, max_episode=500, layer_num=3, layer_dim=200, seed=1234, ctrl0=[], load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      rho0:
     |         --description: initial state (density matrix).
     |         --type: matrix
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |      
     |      savefile:
     |          --description: True: save the control coefficients and the value of the target function
     |                               for each episode.
     |                         False: save the control coefficients and all the value of the target
     |                                function for the last episode.
     |          --type: bool
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  mintime(self, f, W=[], M=[], method='binary', target='QFIM', LDtype='SLD')
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc, decay=[], ctrl_bound=[])
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class DDPG_Sopt(quanestimation.StateOpt.StateStruct.StateSystem)
     |  DDPG_Sopt(savefile=False, max_episode=500, layer_num=3, layer_dim=200, seed=1234, psi0=[], load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      DDPG_Sopt
     |      quanestimation.StateOpt.StateStruct.StateSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, max_episode=500, layer_num=3, layer_dim=200, seed=1234, psi0=[], load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      savefile:
     |          --description: True: save the states and the value of target function for each episode .
     |                         False: save the states and the value of target function for the last episode.
     |          --type: bool
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      seed:
     |          --description: random seed.
     |          --type: int
     |      
     |      eps:
     |          --description: machine eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class DE_Compopt(quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem)
     |  DE_Compopt(savefile=False, popsize=10, psi0=[], ctrl0=[], measurement0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, eps=1e-08)
     |  
     |  Method resolution order:
     |      DE_Compopt
     |      quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CM(self, rho0, W=[])
     |      Description: use DE algorithm to optimize control coefficients and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      rho0:
     |          --description: initial state.
     |          --type: density matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SC(self, W=[], M=[], target='QFIM', LDtype='SLD')
     |      Description: use DE algorithm to optimize states and control coefficients.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SCM(self, W=[])
     |      Description: use DE algorithm to optimize states, control coefficients and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SM(self, W=[])
     |      Description: use DE algorithm to optimize states and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, popsize=10, psi0=[], ctrl0=[], measurement0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      
     |      psi0:
     |         --description: initial guesses of states (kets).
     |         --type: array
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      measurement0:
     |         --description: a set of POVMs.
     |         --type: list (of vector)
     |      
     |      savefile:
     |          --description: True: save the states (or controls, measurements) and the value of the
     |                               target function for each episode.
     |                         False: save the states (or controls, measurements) and all the value
     |                                 of the target function for the last episode.
     |          --type: bool
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[], ctrl_bound=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      psi0:
     |         --description: initial state.
     |         --type: vector
     |      
     |      measurement0:
     |         --description: a set of POVMs.
     |         --type: list (of vector)
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save_meas(self)
     |  
     |  load_save_state(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class DE_Copt(quanestimation.ControlOpt.ControlStruct.ControlSystem)
     |  DE_Copt(savefile=False, popsize=10, ctrl0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      DE_Copt
     |      quanestimation.ControlOpt.ControlStruct.ControlSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   HCRB.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, popsize=10, ctrl0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      rho0:
     |         --description: initial state (density matrix).
     |         --type: matrix
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |      
     |      savefile:
     |          --description: True: save the control coefficients and the value of the target function
     |                               for each episode.
     |                         False: save the control coefficients and all the value of the target
     |                                function for the last episode.
     |          --type: bool
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  mintime(self, f, W=[], M=[], method='binary', target='QFIM', LDtype='SLD')
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc, decay=[], ctrl_bound=[])
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class DE_Mopt(quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem)
     |  DE_Mopt(mtype, minput, savefile=False, popsize=10, measurement0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      DE_Mopt
     |      quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, W=[])
     |      Description: use differential evolution algorithm to update the measurements that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, mtype, minput, savefile=False, popsize=10, measurement0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, load=False, eps=1e-08)
     |       ----------
     |       Inputs
     |       ----------
     |      savefile:
     |           --description: True: save the measurements and the value of the target function
     |                                for each episode.
     |                          False: save the measurements and all the value of the target
     |                                 function for the last episode.
     |           --type: bool
     |      
     |       measurement0:
     |          --description: a set of POVMs.
     |          --type: list (of vector)
     |      
     |       eps:
     |           --description: calculation eps.
     |           --type: float
     |      
     |       notes: the Weyl-Heisenberg covariant SIC-POVM fiducial state of dimension $d$
     |              are download from http://www.physics.umb.edu/Research/QBism/solutions.html.
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, rho0, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class DE_Sopt(quanestimation.StateOpt.StateStruct.StateSystem)
     |  DE_Sopt(savefile=False, popsize=10, psi0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      DE_Sopt
     |      quanestimation.StateOpt.StateStruct.StateSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, popsize=10, psi0=[], max_episode=1000, c=1.0, cr=0.5, seed=1234, load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      savefile:
     |          --description: True: save the states and the value of target function for each episode .
     |                         False: save the states and the value of target function for the last episode.
     |          --type: bool
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      seed:
     |          --description: random seed.
     |          --type: int
     |      
     |      eps:
     |          --description: machine eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class GRAPE_Copt(quanestimation.ControlOpt.ControlStruct.ControlSystem)
     |  GRAPE_Copt(savefile=False, Adam=True, ctrl0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, load=False, eps=1e-08, auto=True)
     |  
     |  Method resolution order:
     |      GRAPE_Copt
     |      quanestimation.ControlOpt.ControlStruct.ControlSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   HCRB.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, Adam=True, ctrl0=[], max_episode=300, epsilon=0.01, beta1=0.9, beta2=0.99, load=False, eps=1e-08, auto=True)
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      rho0:
     |         --description: initial state (density matrix).
     |         --type: matrix
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |      
     |      savefile:
     |          --description: True: save the control coefficients and the value of the target function
     |                               for each episode.
     |                         False: save the control coefficients and all the value of the target
     |                                function for the last episode.
     |          --type: bool
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  mintime(self, f, W=[], M=[], method='binary', target='QFIM', LDtype='SLD')
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc, decay=[], ctrl_bound=[])
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class Lindblad(builtins.object)
     |  Lindblad(tspan, rho0, H0, dH, decay=[], Hc=[], ctrl=[])
     |  
     |  General dynamics of density matrices in the form of time local Lindblad master equation.
     |  {\partial_t 
ho} = -i[H, 
ho] + \sum_n {\gamma_n} {Ln.rho.Ln^{\dagger}
     |               -0.5(rho.Ln^{\dagger}.Ln+Ln^{\dagger}.Ln.rho)}.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, tspan, rho0, H0, dH, decay=[], Hc=[], ctrl=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      rho0:
     |         --description: initial state (density matrix).
     |         --type: matrix
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |         --description: control coefficients.
     |         --type: list (of array)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0] represent a list of decay operators and
     |                        decay[1] represent the corresponding decay rates.
     |         --type: list
     |  
     |  expm(self)
     |  
     |  secondorder_derivative(self, d2H)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class NM_Sopt(quanestimation.StateOpt.StateStruct.StateSystem)
     |  NM_Sopt(savefile=False, state_num=10, psi0=[], max_episode=1000, ar=1.0, ae=2.0, ac=0.5, as0=0.5, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      NM_Sopt
     |      quanestimation.StateOpt.StateStruct.StateSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, state_num=10, psi0=[], max_episode=1000, ar=1.0, ae=2.0, ac=0.5, as0=0.5, seed=1234, load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      savefile:
     |          --description: True: save the states and the value of target function for each episode .
     |                         False: save the states and the value of target function for the last episode.
     |          --type: bool
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      seed:
     |          --description: random seed.
     |          --type: int
     |      
     |      eps:
     |          --description: machine eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class PSO_Compopt(quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem)
     |  PSO_Compopt(savefile=False, particle_num=10, psi0=[], ctrl0=[], measurement0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, eps=1e-08)
     |  
     |  Method resolution order:
     |      PSO_Compopt
     |      quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CM(self, rho0, W=[])
     |      Description: use DE algorithm to optimize control coefficients and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      rho0:
     |          --description: initial state.
     |          --type: density matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SC(self, W=[], M=[], target='QFIM', LDtype='SLD')
     |      Description: use DE algorithm to optimize states and control coefficients.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SCM(self, W=[])
     |      Description: use DE algorithm to optimize states, control coefficients and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  SM(self, W=[])
     |      Description: use DE algorithm to optimize states and the measurements.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, particle_num=10, psi0=[], ctrl0=[], measurement0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      
     |      psi0:
     |         --description: initial guesses of states (kets).
     |         --type: array
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      measurement0:
     |         --description: a set of POVMs.
     |         --type: list (of vector)
     |      
     |      savefile:
     |          --description: True: save the states (or controls, measurements) and the value of the
     |                               target function for each episode.
     |                         False: save the states (or controls, measurements) and all the value
     |                                 of the target function for the last episode.
     |          --type: bool
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[], ctrl_bound=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      psi0:
     |         --description: initial state.
     |         --type: vector
     |      
     |      measurement0:
     |         --description: a set of POVMs.
     |         --type: list (of vector)
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save_meas(self)
     |  
     |  load_save_state(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ComprehensiveOpt.ComprehensiveStruct.ComprehensiveSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class PSO_Copt(quanestimation.ControlOpt.ControlStruct.ControlSystem)
     |  PSO_Copt(savefile=False, particle_num=10, ctrl0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      PSO_Copt
     |      quanestimation.ControlOpt.ControlStruct.ControlSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   HCRB.
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use differential evolution algorithm to update the control coefficients that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, particle_num=10, ctrl0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      rho0:
     |         --description: initial state (density matrix).
     |         --type: matrix
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix or a list of matrix
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |      
     |      savefile:
     |          --description: True: save the control coefficients and the value of the target function
     |                               for each episode.
     |                         False: save the control coefficients and all the value of the target
     |                                function for the last episode.
     |          --type: bool
     |      
     |      ctrl0:
     |          --description: initial control coefficients.
     |          --type: list (of vector)
     |      
     |      eps:
     |          --description: calculation eps.
     |          --type: float
     |  
     |  mintime(self, f, W=[], M=[], method='binary', target='QFIM', LDtype='SLD')
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc, decay=[], ctrl_bound=[])
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.ControlOpt.ControlStruct.ControlSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class PSO_Mopt(quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem)
     |  PSO_Mopt(mtype, minput, savefile=False, particle_num=10, measurement0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      PSO_Mopt
     |      quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, W=[])
     |      Description: use differential evolution algorithm to update the measurements that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, mtype, minput, savefile=False, particle_num=10, measurement0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, load=False, eps=1e-08)
     |       ----------
     |       Inputs
     |       ----------
     |      savefile:
     |           --description: True: save the measurements and the value of the target function
     |                                for each episode.
     |                          False: save the measurements and all the value of the target
     |                                 function for the last episode.
     |           --type: bool
     |      
     |       measurement0:
     |          --description: a set of POVMs.
     |          --type: list (of vector)
     |      
     |       eps:
     |           --description: calculation eps.
     |           --type: float
     |      
     |       notes: the Weyl-Heisenberg covariant SIC-POVM fiducial state of dimension $d$
     |              are download from http://www.physics.umb.edu/Research/QBism/solutions.html.
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem:
     |  
     |  dynamics(self, tspan, rho0, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, rho0, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.MeasurementOpt.MeasurementStruct.MeasurementSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class PSO_Sopt(quanestimation.StateOpt.StateStruct.StateSystem)
     |  PSO_Sopt(savefile=False, particle_num=10, psi0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, load=False, eps=1e-08)
     |  
     |  Method resolution order:
     |      PSO_Sopt
     |      quanestimation.StateOpt.StateStruct.StateSystem
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[])
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   CFI (1/Tr(WF^{-1} with F the CFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      M:
     |          --description: a set of POVM.
     |          --type: list of matrix
     |      
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  HCRB(self, W=[])
     |  
     |  QFIM(self, W=[], LDtype='SLD')
     |      Description: use autodifferential algorithm to search the optimal initial state that maximize the
     |                   QFI (1/Tr(WF^{-1} with F the QFIM).
     |      
     |      ---------
     |      Inputs
     |      ---------
     |      W:
     |          --description: weight matrix.
     |          --type: matrix
     |  
     |  __init__(self, savefile=False, particle_num=10, psi0=[], max_episode=[1000, 100], c0=1.0, c1=2.0, c2=2.0, seed=1234, load=False, eps=1e-08)
     |      ----------
     |      Inputs
     |      ----------
     |      savefile:
     |          --description: True: save the states and the value of target function for each episode .
     |                         False: save the states and the value of target function for the last episode.
     |          --type: bool
     |      
     |      psi0:
     |          --description: initial guess of states (kets).
     |          --type: array
     |      
     |      seed:
     |          --description: random seed.
     |          --type: int
     |      
     |      eps:
     |          --description: machine eps.
     |          --type: float
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  dynamics(self, tspan, H0, dH, Hc=[], ctrl=[], decay=[])
     |      ----------
     |      Inputs
     |      ----------
     |      tspan:
     |         --description: time series.
     |         --type: array
     |      
     |      H0:
     |         --description: free Hamiltonian.
     |         --type: matrix (a list of matrix)
     |      
     |      dH:
     |         --description: derivatives of Hamiltonian on all parameters to
     |                        be estimated. For example, dH[0] is the derivative
     |                        vector on the first parameter.
     |         --type: list (of matrix)
     |      
     |      Hc:
     |         --description: control Hamiltonian.
     |         --type: list (of matrix)
     |      
     |      ctrl:
     |          --description: control coefficients.
     |          --type: list (of vector)
     |      
     |      decay:
     |         --description: decay operators and the corresponding decay rates.
     |                        decay[0][0] represent the first decay operator and
     |                        decay[0][1] represent the corresponding decay rate.
     |         --type: list
     |      
     |      ctrl_bound:
     |         --description: lower and upper bounds of the control coefficients.
     |                        ctrl_bound[0] represent the lower bound of the control coefficients and
     |                        ctrl_bound[1] represent the upper bound of the control coefficients.
     |         --type: list
     |  
     |  kraus(self, K, dK)
     |  
     |  load_save(self)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from quanestimation.StateOpt.StateStruct.StateSystem:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class adaptMZI(builtins.object)
     |  adaptMZI(x, p, rho0)
     |  
     |  Methods defined here:
     |  
     |  __init__(self, x, p, rho0)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  general(self)
     |  
     |  offline(self, method='DE', popsize=10, particle_num=10, DeltaPhi0=[], c=1.0, cr=0.5, c0=1.0, c1=2.0, c2=2.0, seed=1234, max_episode=1000, eps=1e-08)
     |  
     |  online(self, output='phi')
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class adaptive(builtins.object)
     |  adaptive(x, p, rho0, max_episode=1000, eps=1e-08)
     |  
     |  Methods defined here:
     |  
     |  CFIM(self, M=[], W=[], savefile=False)
     |  
     |  Mopt(self, W=[])
     |  
     |  __init__(self, x, p, rho0, max_episode=1000, eps=1e-08)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  dynamics(self, tspan, H, dH, Hc=[], ctrl=[], decay=[])
     |  
     |  kraus(self, K, dK)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)

FUNCTIONS
    AdaptiveInput(x, func, dfunc, channel='dynamics')
    
    BCFIM(x, p, rho, drho, M=[], eps=1e-08)
        Description: Calculation Bayesian version of classical Fisher information
                     matrix (CFIM) for a density matrix.
        
        ---------
        Inputs
        ---------
        x:
            --description: the regimes of x for the integral.
            --type: list of arrays
        
        p:
            --description: the prior distribution.
            --type: multidimensional array
        
        rho:
            --description: parameterized density matrix.
            --type: multidimensional lists
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated.
            --type: multidimensional lists
        
        M:
           --description: a set of POVM. It takes the form [M1, M2, ...].
           --type: list (of matrix)
        
        ----------
        Returns
        ----------
        BCFIM:
            --description: Bayesian version of classical Fisher information matrix.
                           If the length of x is one, the output is a float number (BCFI),
                           otherwise it returns a matrix (BCFIM).
            --type: float number (BCFI) or matrix (BCFIM)
    
    BCRB(x, p, rho, drho, M=[], b=[], db=[], btype=1, eps=1e-08)
    
    BQCRB(x, p, rho, drho, b=[], db=[], btype=1, LDtype='SLD', eps=1e-08)
    
    BQFIM(x, p, rho, drho, LDtype='SLD', eps=1e-08)
        Description: Calculation of Bayesian version of quantum Fisher information
                     matrix (QFIM) for a density matrix.
        
        ----------
        Inputs
        ----------
        x:
            --description: the regimes of x for the integral.
            --type: list of arrays
        
        p:
            --description: the prior distribution.
            --type: multidimensional array
        
        rho:
            --description: parameterized density matrix.
            --type: multidimensional lists
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example,
            --type: multidimensional lists
        
        LDtype:
            --description: the type of logarithmic derivatives.
            --type: string {'SLD', 'RLD', 'LLD'}
        
        ----------
        Returns
        ----------
        BQFIM:
            --description: Bayesian version of quantum Fisher information matrix.
                           If the length of x is 1, the output is a float number (BQFI),
                           otherwise it returns a matrix (BQFIM).
            --type: float number (BQFI) or matrix (BQFIM)
    
    Bayes(x, p, rho, y, M=[], savefile=False)
    
    CFIM(rho, drho, M=[], eps=1e-08)
        Description: Calculation classical Fisher information matrix (CFIM)
                     for a density matrix.
        
        ---------
        Inputs
        ---------
        rho:
            --description: parameterized density matrix.
            --type: matrix
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example, drho[0] is the derivative
                           vector on the first parameter.
            --type: list (of matrix)
        
        M:
           --description: a set of POVM. It takes the form [M1, M2, ...].
           --type: list (of matrix)
        
        ----------
        Returns
        ----------
        CFIM:
            --description: classical Fisher information matrix. If the length
                           of drho is one, the output is a float number (CFI),
                           otherwise it returns a matrix (CFIM).
            --type: float number (CFI) or matrix (CFIM)
    
    ComprehensiveOpt(savefile=False, method='AD', **kwargs)
    
    ControlOpt(savefile=False, method='auto-GRAPE', **kwargs)
    
    FIM(p, dp, eps=1e-08)
    
    HCRB(rho, drho, W, eps=1e-08)
        Description: solve Holevo Cramer-Rao bound via the semidefinite program (SDP).
        
        ---------
        Inputs
        ---------
        rho:
            --description: parameterized density matrix.
            --type: matrix
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example, drho[0] is the derivative
                           vector on the first parameter.
            --type: list (of matrix)
        
        W:
           --description: cost matrix.
           --type: matrix
        
        ----------
        Returns
        ----------
        Holevo bound:
            --description: the Holevo Cramer-Rao bound.
            --type: float
    
    LLD(rho, drho, rep='original', eps=1e-08)
        Description: Calculation of the left logarithmic derivative (LLD)
                    for a density matrix.
        
        ----------
        Inputs
        ----------
        rho:
            --description: parameterized density matrix.
            --type: matrix
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example, drho[0] is the derivative
                           vector on the first parameter.
            --type: list (of matrix)
        
        rep:
            --description: the basis for the LLDs.
                           rep=original means the basis for obtained LLDs is the
                           same with the density matrix (rho).
                           rep=eigen means the LLDs are written in the eigenspace of
                           the density matrix (rho).
            --type: string {"original", "eigen"}
        
        ----------
        Returns
        ----------
        LLD:
            --description: LLD for the density matrix (rho).
            --type: list (of matrix)
    
    MLE(x, rho, y, M=[], savefile=False)
    
    MeasurementOpt(mtype='projection', minput=[], savefile=False, method='DE', **kwargs)
    
    OBB(x, p, dp, rho, drho, d2rho, LDtype='SLD', eps=1e-08)
    
    QFIM(rho, drho, LDtype='SLD', exportLD=False, eps=1e-08)
        Description: Calculation of quantum Fisher information matrix (QFIM)
                    for a density matrix.
        
        ----------
        Inputs
        ----------
        rho:
            --description: parameterized density matrix.
            --type: matrix
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example, drho[0] is the derivative
                           vector on the first parameter.
            --type: list (of matrix)
        
        LDtype:
            --description: the type of logarithmic derivatives.
            --type: string {'SLD', 'RLD', 'LLD'}
        
        exportLD:
            --description: if True, the corresponding value of logarithmic derivatives
                           will be exported.
            --type: bool
        
        ----------
        Returns
        ----------
        QFIM:
            --description: Quantum Fisher information matrix. If the length
                           of drho is 1, the output is a float number (QFI),
                           otherwise it returns a matrix (QFIM).
            --type: float number (QFI) or matrix (QFIM)
    
    QFIM_Bloch(r, dr, eps=1e-08)
        Description: Calculation of quantum Fisher information matrix (QFIM)
                    in Bloch representation.
        
        ----------
        Inputs
        ----------
        r:
            --description: parameterized Bloch vector.
            --type: vector
        
        dr:
            --description: derivatives of Bloch vector on all parameters to
                            be estimated. For example, dr[0] is the derivative
                            vector on the first parameter.
            --type: list (of vector)
    
    QFIM_Gauss(R, dR, D, dD, eps=1e-08)
        Description: Calculation of quantum Fisher information matrix (QFIM)
                    for Gaussian states representation.
        
        ----------
        Inputs
        ----------
        R:
        
        dR:
        
        D:
        
        dD:
    
    QFIM_Kraus(rho0, K, dK, eps=1e-08)
    
    QVTB(x, p, dp, rho, drho, btype=1, LDtype='SLD', eps=1e-08)
    
    QZZB(x, p, rho, eps=1e-08)
    
    RLD(rho, drho, rep='original', eps=1e-08)
        Description: calculation of the right logarithmic derivative (RLD)
                     for a density matrix.
        
        ----------
        Inputs
        ----------
        rho:
            --description: parameterized density matrix.
            --type: matrix
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example, drho[0] is the derivative
                           vector on the first parameter.
            --type: list (of matrix)
        
        rep:
            --description: the basis for the RLDs.
                           rep=original means the basis for obtained RLDs is the
                           same with the density matrix (rho).
                           rep=eigen means the RLDs are written in the eigenspace of
                           the density matrix (rho).
            --type: string {"original", "eigen"}
        
        ----------
        Returns
        ----------
        RLD:
            --description: RLD for the density matrix (rho).
            --type: list (of matrix)
    
    SIC(dim)
    
    SLD(rho, drho, rep='original', eps=1e-08)
        Description: calculation of the symmetric logarithmic derivative (SLD)
                     for a density matrix.
        
        ----------
        Inputs
        ----------
        rho:
            --description: parameterized density matrix.
            --type: matrix
        
        drho:
            --description: derivatives of density matrix (rho) on all parameters
                           to be estimated. For example, drho[0] is the derivative
                           vector on the first parameter.
            --type: list (of matrix)
        
        rep:
            --description: the basis for the SLDs.
                           rep=original means the basis for obtained SLDs is the
                           same with the density matrix (rho).
                           rep=eigen means the SLDs are written in the eigenspace of
                           the density matrix (rho).
            --type: string {"original", "eigen"}
        
        ----------
        Returns
        ----------
        SLD:
            --description: SLD for the density matrix (rho).
            --type: list (of matrix)
    
    SpinSqueezing(rho, basis='Dicke', output='KU')
        Description: Calculate spin squeezing parameter for a density matrix.
        
        ---------
        Inputs
        ---------
        N :
           --description: particle number.
           --type: int
        rho :
           --description: density matrix.
           --type: matrix
        
        output:
           --description: if output=='KU',the output of the squeezing parameter is defined by
                          Kitagawa and Ueda. if output=='WBIMH',the output of the squeezing
                          parameter is defined by Wineland et al.
           --type: string {'KU', 'WBIMH'}
        
        ----------
        Returns
        ----------
        Xi :
           --description: squeezing_parameter.
           --type: float
    
    StateOpt(savefile=False, method='AD', **kwargs)
    
    TargetTime(f, tspan, func, *args, **kwargs)
    
    VTB(x, p, dp, rho, drho, M=[], btype=1, eps=1e-08)
    
    annihilation(n)
    
    basis(dim, index)
    
    csv2npy_controls(controls, num)
    
    csv2npy_measurements(M, num)
    
    csv2npy_states(states, num=1)
    
    gramschmidt(A)
    
    kraus(rho0, K, dK)
    
    mat_vec_convert(A)
    
    suN_generator(n)

DATA
    __all__ = ['ControlOpt', 'StateOpt', 'MeasurementOpt', 'ComprehensiveO...

FILE
    /workspaces/QuanEstimation/quanestimation/__init__.py


